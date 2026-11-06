
#include <Rcpp.h>
#include <thread>
#include <vector>
#include <cmath>
#include <algorithm>
using namespace Rcpp;

inline double gh_l0_clamp(double x, double tiny) {
  if (!std::isfinite(x) || x < tiny) return tiny;
  return x;
}

// [[Rcpp::export]]
Rcpp::List gh_L0_ZIP_block_cpp(
    Rcpp::List PNy_list,
    Rcpp::List PN0_list,
    Rcpp::List dPNy_list,
    Rcpp::List dPN0_list,
    Rcpp::List pb_list,
    Rcpp::List y0_list,
    Rcpp::LogicalVector zero_inflation,
    Rcpp::NumericVector WN,
    int J,
    bool need_grad,
    double tiny = 1e-300,
    int n_threads = 1) {

  const int k = PNy_list.size();
  if (k < 1) stop("Empty PNy_list.");
  NumericMatrix P0 = PNy_list[0];
  const int nb = P0.nrow();
  const int QN = P0.ncol();
  if (WN.size() != QN) stop("WN length does not match GH node count.");
  if (zero_inflation.size() != k) stop("zero_inflation length mismatch.");

  std::vector<const double*> Py(k), Pz(k, nullptr), PB(k, nullptr);
  std::vector<const int*> Y0(k, nullptr);
  std::vector<const double*> dPy(k, nullptr), dPz(k, nullptr);

  for (int j=0; j<k; ++j) {
    SEXP sy = PNy_list[j];
    NumericMatrix my(sy);
    if (my.nrow()!=nb || my.ncol()!=QN) stop("PNy dimension mismatch.");
    Py[j] = REAL(sy);
    if (zero_inflation[j]) {
      SEXP sz = PN0_list[j];
      NumericMatrix mz(sz);
      if (mz.nrow()!=nb || mz.ncol()!=QN) stop("PN0 dimension mismatch.");
      Pz[j] = REAL(sz);
      SEXP sp = pb_list[j];
      NumericVector vp(sp);
      if (vp.size()!=nb) stop("pb length mismatch.");
      PB[j] = REAL(sp);
      SEXP s0 = y0_list[j];
      LogicalVector v0(s0);
      if (v0.size()!=nb) stop("y0 length mismatch.");
      Y0[j] = LOGICAL(s0);
    }
    if (need_grad) {
      SEXP sd = dPNy_list[j];
      NumericMatrix md(sd);
      if (md.nrow()!=nb || md.ncol()!=QN*J) stop("dPNy dimension mismatch.");
      dPy[j] = REAL(sd);
      if (zero_inflation[j]) {
        SEXP sd0 = dPN0_list[j];
        NumericMatrix md0(sd0);
        if (md0.nrow()!=nb || md0.ncol()!=QN*J) stop("dPN0 dimension mismatch.");
        dPz[j] = REAL(sd0);
      }
    }
  }

  int nt = std::max(1, std::min(n_threads, nb));
  std::vector<double> nll_thr(nt, 0.0);
  std::vector< std::vector<double> > grad_thr;
  if (need_grad) grad_thr.assign(nt, std::vector<double>(k*J, 0.0));
  std::vector<int> bad_thr(nt, 0);

  auto worker = [&](int begin, int end, int tid) {
    std::vector<double> T(k), prefix(k+1), suffix(k+1);
    std::vector<double> dL;
    if (need_grad) dL.resize(k*J);

    for (int i=begin; i<end; ++i) {
      double Li = 0.0;
      if (need_grad) std::fill(dL.begin(), dL.end(), 0.0);

      for (int q=0; q<QN; ++q) {
        for (int j=0; j<k; ++j) {
          const double py = gh_l0_clamp(Py[j][i + nb*q], tiny);
          double tj;
          if (!zero_inflation[j]) {
            tj = py;
          } else {
            const double pb = PB[j][i];
            if (Y0[j][i]) {
              const double pz = gh_l0_clamp(Pz[j][i + nb*q], tiny);
              tj = (1.0 - pb) + pb * pz;
            } else {
              tj = pb * py;
            }
          }
          T[j] = gh_l0_clamp(tj, tiny);
        }

        prefix[0] = 1.0;
        for (int j=0; j<k; ++j) prefix[j+1] = prefix[j] * T[j];
        suffix[k] = 1.0;
        for (int j=k-1; j>=0; --j) suffix[j] = suffix[j+1] * T[j];

        const double w = WN[q];
        const double prod = prefix[k];
        Li += w * prod;

        if (need_grad) {
          for (int j=0; j<k; ++j) {
            const double except = prefix[j] * suffix[j+1];
            const double pb = zero_inflation[j] ? PB[j][i] : 1.0;
            for (int r=0; r<J; ++r) {
              const int col = r*QN + q;
              double dt;
              if (!zero_inflation[j]) {
                dt = dPy[j][i + nb*col];
              } else if (Y0[j][i]) {
                dt = pb * dPz[j][i + nb*col];
              } else {
                dt = pb * dPy[j][i + nb*col];
              }
              dL[j*J+r] += w * except * dt;
            }
          }
        }
      }

      if (!std::isfinite(Li) || Li <= tiny) {
        bad_thr[tid] = 1;
        continue;
      }
      nll_thr[tid] -= std::log(Li);
      if (need_grad) {
        for (int g=0; g<k*J; ++g) grad_thr[tid][g] -= dL[g] / Li;
      }
    }
  };

  std::vector<std::thread> threads;
  threads.reserve(nt);
  const int chunk = (nb + nt - 1) / nt;
  for (int t=0; t<nt; ++t) {
    const int begin = t*chunk;
    const int end = std::min(nb, begin+chunk);
    if (begin >= end) break;
    threads.emplace_back(worker, begin, end, t);
  }
  for (auto &th : threads) th.join();

  for (int t=0; t<nt; ++t) if (bad_thr[t]) {
    return List::create(_["nll"] = R_PosInf,
                        _["grad_Psi"] = R_NilValue,
                        _["threads_used"] = nt);
  }

  double nll = 0.0;
  for (int t=0; t<nt; ++t) nll += nll_thr[t];
  NumericMatrix grad(k, J);
  if (need_grad) {
    for (int t=0; t<nt; ++t)
      for (int j=0; j<k; ++j)
        for (int r=0; r<J; ++r)
          grad(j,r) += grad_thr[t][j*J+r];
  }

  return List::create(_["nll"] = nll,
                      _["grad_Psi"] = need_grad ? (SEXP)grad : R_NilValue,
                      _["threads_used"] = nt);
}
