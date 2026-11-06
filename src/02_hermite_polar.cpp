#include <RcppArmadillo.h>
#include <algorithm>
#include <limits>
#include <vector>
#include <thread>
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

namespace zipm_parallel {
inline int effective_threads(int requested, int nwork) {
    if (nwork <= 1) return 1;
    int nt = std::max(1, requested);
    nt = std::min(nt, nwork);
    return nt;
}

template <class Func>
inline int parallel_chunks(int nwork, int requested, Func fn) {
    const int nt = effective_threads(requested, nwork);
    if (nt <= 1) { fn(0, nwork, 0); return 1; }
    const int chunk = (nwork + nt - 1) / nt;
    std::vector<std::thread> workers;
    workers.reserve(nt);
    int launched = 0;
    for (int t = 0; t < nt; ++t) {
        const int begin = t * chunk;
        const int end = std::min(nwork, begin + chunk);
        if (begin >= end) break;
        workers.emplace_back([&, begin, end, t]() { fn(begin, end, t); });
        ++launched;
    }
    for (auto& th : workers) th.join();
    return launched;
}
}

// -----------------------------------------------------------------------------
// Hermite polynomial kernels
//  * never constructs k x k diagonal matrices
//  * uses one Cholesky per observation (no second SPD factorization)
//  * performs normalized fixed-grid reweighting entirely in C++
//  * does not materialize Abar / score matrices
// -----------------------------------------------------------------------------

// [[Rcpp::export]]
Rcpp::List hermite_state_cpp(
    const arma::mat& Psi,   // k x J
    const arma::mat& C0,    // nb x k
    const arma::mat& C1,    // nb x k
    const arma::mat& C2,    // nb x k
    const double tiny_r = 1e-10,
    const int n_threads = 1
) {
    const int k = Psi.n_rows;
    const int J = Psi.n_cols;
    const int nb = C0.n_rows;

    arma::vec r(k);
    arma::mat U(k, J, arma::fill::zeros);
    for (int j = 0; j < k; ++j) {
        const double rr2 = arma::dot(Psi.row(j), Psi.row(j));
        r(j) = std::sqrt(rr2);
        if (!std::isfinite(r(j)) || r(j) <= tiny_r) {
            Rcpp::stop("hermite_state_cpp: loading-row norm <= tiny_r.");
        }
        U.row(j) = Psi.row(j) / r(j);
    }

    arma::vec logbase(nb); logbase.fill(NA_REAL);
    arma::mat MUX(nb, k); MUX.fill(NA_REAL);
    arma::mat VARX(nb, k); VARX.fill(NA_REAL);
    arma::mat MU(nb, J); MU.fill(NA_REAL);
    arma::mat VPACK(nb, J * J); VPACK.fill(NA_REAL);

    const arma::mat IJ = arma::eye(J, J);

    const int nt = zipm_parallel::effective_threads(n_threads, nb);

    // Each observation is independent. std::thread is used deliberately so
    // n_threads works even when R/Rtools is built without OpenMP.
    zipm_parallel::parallel_chunks(nb, nt, [&](int i_begin, int i_end, int /*tid*/) {
    for (int ii = i_begin; ii < i_end; ++ii) {
        arma::vec l(k), q(k);
        double cst = 0.0;
        bool valid = true;

        for (int j = 0; j < k; ++j) {
            const double rr = r(j);
            const double c0 = C0(ii, j);
            const double c1 = C1(ii, j);
            const double c2 = C2(ii, j);
            l(j) = c1 / rr;
            q(j) = c2 / (rr * rr);
            cst += c0 - c2;
            if (!std::isfinite(l(j)) || !std::isfinite(q(j)) || !std::isfinite(c0)) {
                valid = false;
                break;
            }
        }
        if (!valid) continue;

        const arma::vec a = Psi.t() * l;

        // B = 2 Psi' diag(q) Psi without forming either diag(q) or k x k data.
        arma::mat weightedPsi = Psi;
        weightedPsi.each_col() %= q;
        arma::mat B = 2.0 * (Psi.t() * weightedPsi);
        arma::mat A = IJ - 0.5 * (B + B.t());

        arma::mat R;
        if (!arma::chol(R, A)) continue;

        // Reuse Cholesky for the inverse; no second SPD factorization.
        arma::mat Rinv;
        if (!arma::solve(Rinv, arma::trimatu(R), IJ, arma::solve_opts::fast)) continue;
        if (!Rinv.is_finite()) continue;
        arma::mat V = Rinv * Rinv.t();

        const arma::vec mu = V * a;
        const double logdetA = 2.0 * arma::sum(arma::log(R.diag()));
        const double lb = cst - 0.5 * logdetA + 0.5 * arma::dot(a, mu);
        if (!std::isfinite(lb)) continue;

        logbase(ii) = lb;
        MU.row(ii) = mu.t();

        int cc = 0;
        for (int rr = 0; rr < J; ++rr)
            for (int ss = 0; ss < J; ++ss)
                VPACK(ii, cc++) = V(rr, ss);

        const arma::vec mux = U * mu;
        const arma::mat UV = U * V;
        const arma::vec varx = arma::sum(UV % U, 1);
        MUX.row(ii) = mux.t();
        VARX.row(ii) = varx.t();
    }
    });

    return Rcpp::List::create(
        Rcpp::Named("logbase") = logbase,
        Rcpp::Named("mux") = MUX,
        Rcpp::Named("varx") = VARX,
        Rcpp::Named("mu") = MU,
        Rcpp::Named("Vpack") = VPACK
    );
}


// [[Rcpp::export]]
Rcpp::List hermite_reweight_cpp(
    const Rcpp::List& LH_list,
    const Rcpp::List& dLH_list,
    const Rcpp::List& FPB_list,
    const arma::mat& MUX,
    const arma::mat& VARX,
    const arma::vec& logbase,
    const arma::mat& C0,
    const arma::mat& C1,
    const arma::mat& C2,
    const arma::mat& D0,
    const arma::mat& D1,
    const arma::mat& D2,
    const arma::vec& herm_x,
    const arma::vec& herm_w,
    const bool need_grad = false,
    const bool estimate_B_factor = false,
    const double tiny_v = 1e-14,
    const int n_threads = 1
) {
    const int nb = MUX.n_rows;
    const int k = MUX.n_cols;
    const int Q = herm_x.n_elem;
    if ((int)LH_list.size() != k) Rcpp::stop("hermite_reweight_cpp: LH_list length != k.");
    if ((int)herm_w.n_elem != Q) Rcpp::stop("hermite_reweight_cpp: herm_x/herm_w size mismatch.");

    arma::mat CORR(nb,k,arma::fill::zeros), CORR_MU(nb,k,arma::fill::zeros), CORR_VAR(nb,k,arma::fill::zeros), CORR_R(nb,k,arma::fill::zeros);
    arma::mat GPB_DIRECT;
    if (need_grad && estimate_B_factor) GPB_DIRECT.zeros(nb,k);

    arma::vec logw_base(Q);
    for(int q=0;q<Q;++q){
        const double wq=std::max(herm_w(q),std::numeric_limits<double>::min());
        const double x=herm_x(q);
        logw_base(q)=std::log(wq)+0.5*x*x;
    }

    // Extract R handles before the OpenMP region.
    std::vector<Rcpp::NumericMatrix> LHm; LHm.reserve(k);
    std::vector<Rcpp::NumericMatrix> dLHm; dLHm.reserve(k);
    std::vector<Rcpp::NumericMatrix> FPBm; FPBm.reserve(k);
    std::vector<int> hasD(k,0), hasP(k,0);
    for(int j=0;j<k;++j){
        Rcpp::NumericMatrix L=LH_list[j];
        if(L.nrow()!=nb||L.ncol()!=Q) Rcpp::stop("hermite_reweight_cpp: LH dimension mismatch.");
        LHm.push_back(L);
        if(need_grad){
            if(j>=dLH_list.size()||Rf_isNull(dLH_list[j])) Rcpp::stop("hermite_reweight_cpp: missing dLH.");
            Rcpp::NumericMatrix D=dLH_list[j]; dLHm.push_back(D); hasD[j]=1;
        } else dLHm.push_back(Rcpp::NumericMatrix(0,0));
        if(need_grad&&estimate_B_factor&&j<FPB_list.size()&&!Rf_isNull(FPB_list[j])){
            Rcpp::NumericMatrix P=FPB_list[j]; FPBm.push_back(P); hasP[j]=(P.nrow()==nb&&P.ncol()==Q);
        } else FPBm.push_back(Rcpp::NumericMatrix(0,0));
    }

    double max_abs_logmass=0.0,sum_abs_logmass=0.0,max_mass_err=0.0,sum_mass_err=0.0;
    long long ndiag=0;
    const int nt=zipm_parallel::effective_threads(n_threads,k);
    std::vector<double> t_max_log(nt,0.0), t_sum_log(nt,0.0), t_max_err(nt,0.0), t_sum_err(nt,0.0);
    std::vector<long long> t_ndiag(nt,0);
    zipm_parallel::parallel_chunks(k, nt, [&](int j_begin, int j_end, int tid) {
        double local_max_log=0.0, local_sum_log=0.0, local_max_err=0.0, local_sum_err=0.0;
        long long local_ndiag=0;
    for(int j=j_begin;j<j_end;++j){
        const Rcpp::NumericMatrix& LH=LHm[j];
        const Rcpp::NumericMatrix& dLH=dLHm[j];
        const Rcpp::NumericMatrix& FPB=FPBm[j];
        for(int i=0;i<nb;++i){
            if(!std::isfinite(logbase(i))) continue;
            const double mu=MUX(i,j),v=VARX(i,j);
            if(!std::isfinite(mu)||!std::isfinite(v)||v<=tiny_v) continue;
            const double invv=1.0/v,inv2v2=0.5*invv*invv,neg_half_logv=-0.5*std::log(v);
            double amax=-std::numeric_limits<double>::infinity();
            for(int q=0;q<Q;++q){const double x=herm_x(q),dx=x-mu,la=logw_base(q)+neg_half_logv-0.5*dx*dx*invv;if(la>amax)amax=la;}
            if(!std::isfinite(amax)) continue;
            double den=0,ef=0,es_mu=0,efs_mu=0,es_var=0,efs_var=0,edr=0,epb=0;
            for(int q=0;q<Q;++q){
                const double x=herm_x(q),dx=x-mu,la=logw_base(q)+neg_half_logv-0.5*dx*dx*invv,aa=std::exp(la-amax),f=LH(i,q);
                if(!std::isfinite(aa)||!std::isfinite(f)) continue;
                den+=aa;ef+=aa*f;
                if(need_grad){const double sm=dx*invv,sv=(dx*dx-v)*inv2v2;es_mu+=aa*sm;efs_mu+=aa*f*sm;es_var+=aa*sv;efs_var+=aa*f*sv;if(hasD[j])edr+=aa*dLH(i,q);if(hasP[j])epb+=aa*FPB(i,q);}
            }
            if(!std::isfinite(den)||den<=0) continue;
            const double G=ef/den,h2=mu*mu+v-1.0;
            CORR(i,j)=G-C0(i,j)-C1(i,j)*mu-C2(i,j)*h2;
            if(need_grad){
                const double Gmu=efs_mu/den-G*(es_mu/den),Gvar=efs_var/den-G*(es_var/den),Gr=edr/den;
                CORR_MU(i,j)=Gmu-C1(i,j)-2.0*C2(i,j)*mu;
                CORR_VAR(i,j)=Gvar-C2(i,j);
                CORR_R(i,j)=Gr-D0(i,j)-D1(i,j)*mu-D2(i,j)*h2;
                if(hasP[j])GPB_DIRECT(i,j)=epb/den;
            }
            const double lm=amax+std::log(den);
            if(std::isfinite(lm)){
                const double alm=std::abs(lm); local_max_log=std::max(local_max_log,alm); local_sum_log+=alm;
                double me=std::numeric_limits<double>::infinity(); if(std::abs(lm)<50.0) me=std::abs(std::expm1(lm));
                local_max_err=std::max(local_max_err,me); local_sum_err+=me; ++local_ndiag;
            }
        }
    }
        t_max_log[tid]=local_max_log; t_sum_log[tid]=local_sum_log;
        t_max_err[tid]=local_max_err; t_sum_err[tid]=local_sum_err; t_ndiag[tid]=local_ndiag;
    });
    for(int t=0;t<nt;++t){
        max_abs_logmass=std::max(max_abs_logmass,t_max_log[t]); sum_abs_logmass+=t_sum_log[t];
        max_mass_err=std::max(max_mass_err,t_max_err[t]); sum_mass_err+=t_sum_err[t]; ndiag+=t_ndiag[t];
    }
    const double mean_abs_logmass=ndiag>0?sum_abs_logmass/(double)ndiag:R_PosInf;
    const double mean_mass_err=ndiag>0?sum_mass_err/(double)ndiag:R_PosInf;
    return Rcpp::List::create(
        Rcpp::Named("CORR")=CORR,Rcpp::Named("CORR_MU")=CORR_MU,Rcpp::Named("CORR_VAR")=CORR_VAR,Rcpp::Named("CORR_R")=CORR_R,Rcpp::Named("GPB_DIRECT")=GPB_DIRECT,
        Rcpp::Named("h2_logmass_max_abs")=max_abs_logmass,Rcpp::Named("h2_logmass_mean_abs")=mean_abs_logmass,Rcpp::Named("h2_mass_max_error")=max_mass_err,Rcpp::Named("h2_mass_mean_error")=mean_mass_err
    );
}

// [[Rcpp::export]]
Rcpp::List hermite_finish_cpp(
    const arma::mat& Psi,
    const arma::mat& C0,
    const arma::mat& C1,
    const arma::mat& C2,
    const arma::mat& D0,
    const arma::mat& D1,
    const arma::mat& D2,
    const arma::vec& logbase,
    const arma::mat& MUX,
    const arma::mat& VARX,
    const arma::mat& MU,
    const arma::mat& VPACK,
    const arma::mat& CORR,
    const arma::mat& CORR_MU,
    const arma::mat& CORR_VAR,
    const arma::mat& CORR_R,
    const arma::vec& obs_weight,
    const bool need_grad = true,
    const double tiny_r = 1e-10,
    const int n_threads = 1
) {
    const int k = Psi.n_rows;
    const int J = Psi.n_cols;
    const int nb = C0.n_rows;

    arma::vec r(k);
    arma::mat U(k, J, arma::fill::zeros);
    for (int j = 0; j < k; ++j) {
        r(j) = std::sqrt(arma::dot(Psi.row(j), Psi.row(j)));
        if (!std::isfinite(r(j)) || r(j) <= tiny_r)
            Rcpp::stop("hermite_finish_cpp: loading-row norm <= tiny_r.");
        U.row(j) = Psi.row(j) / r(j);
    }

    arma::vec logv = logbase + arma::sum(CORR, 1);
    arma::mat G(k, J, arma::fill::zeros);
    arma::mat adjC1(nb, k, arma::fill::zeros);
    arma::mat adjC2(nb, k, arma::fill::zeros);

    if (!need_grad) {
        return Rcpp::List::create(
            Rcpp::Named("log_value") = logv,
            Rcpp::Named("grad_Psi") = G,
            Rcpp::Named("adj_c1") = adjC1,
            Rcpp::Named("adj_c2") = adjC2
        );
    }

    const int nt = zipm_parallel::effective_threads(n_threads, nb);
    std::vector<arma::mat> Gthread;
    Gthread.reserve(nt);
    for (int t = 0; t < nt; ++t)
        Gthread.emplace_back(k, J, arma::fill::zeros);

    zipm_parallel::parallel_chunks(nb, nt, [&](int i_begin, int i_end, int tid) {
    arma::mat& Glocal = Gthread[tid];
    for (int ii = i_begin; ii < i_end; ++ii) {
        if (!std::isfinite(logbase(ii))) continue;

        const arma::vec mu = MU.row(ii).t();
        arma::mat V(J, J);
        int cc = 0;
        bool ok = true;
        for (int rr = 0; rr < J; ++rr) {
            for (int ss = 0; ss < J; ++ss) {
                const double vv = VPACK(ii, cc++);
                if (!std::isfinite(vv)) ok = false;
                V(rr, ss) = vv;
            }
        }
        if (!ok) continue;

        arma::vec l(k), q(k), dl(k), dq(k);
        for (int j = 0; j < k; ++j) {
            const double rr = r(j);
            const double c1 = C1(ii, j);
            const double c2 = C2(ii, j);
            const double d1 = D1(ii, j);
            const double d2 = D2(ii, j);
            l(j) = c1 / rr;
            q(j) = c2 / (rr * rr);
            dl(j) = d1 / rr - c1 / (rr * rr);
            dq(j) = d2 / (rr * rr) - 2.0 * c2 / (rr * rr * rr);
        }

        const arma::vec fmu = CORR_MU.row(ii).t();
        const arma::vec fvar = CORR_VAR.row(ii).t();
        const arma::vec g_mu = U.t() * fmu;

        // G_V = U' diag(fvar) U without a k x k diagonal matrix.
        arma::mat weightedU = U;
        weightedU.each_col() %= fvar;
        const arma::mat G_V = U.t() * weightedU;

        const arma::vec h = V * g_mu;
        const arma::vec g_a = mu + h;
        const arma::mat G_B =
            0.5 * (V + mu * mu.t()) +
            0.5 * (mu * h.t() + h * mu.t()) +
            V * G_V * V;

        // Two dense k x J by J x J GEMMs replace k separate J x J matvecs.
        const arma::mat PsiGB = Psi * G_B;
        const arma::mat PsiV  = Psi * V;

        arma::mat Gi = l * g_a.t();
        arma::mat qPsiGB = PsiGB;
        qPsiGB.each_col() %= q;
        Gi += 4.0 * qPsiGB;

        const arma::vec psiga = Psi * g_a;

        for (int j = 0; j < k; ++j) {
            const double rr = r(j);
            const double invr = 1.0 / rr;
            const arma::rowvec uj = U.row(j);

            // psi' G_B psi and u' G_B u from the already-computed PsiGB.
            const double psiGBpsi = arma::dot(PsiGB.row(j), Psi.row(j));
            const double uGBu = psiGBpsi * invr * invr;

            const double gr =
                D0(ii, j) - D2(ii, j) +
                psiga(j) * dl(j) +
                2.0 * psiGBpsi * dq(j) +
                CORR_R(ii, j);

            // V u_j = V psi_j / r_j; use the PsiV GEMM rather than a matvec.
            arma::vec gu = fmu(j) * mu +
                (2.0 * fvar(j) * invr) * PsiV.row(j).t();
            const double udotgu = arma::dot(uj, gu);
            const arma::vec gu_tangent = (gu - uj.t() * udotgu) * invr;
            Gi.row(j) += gr * uj + gu_tangent.t();

            const double mux = MUX(ii, j);
            const double varx = VARX(ii, j);
            adjC1(ii, j) = psiga(j) * invr - mux;
            adjC2(ii, j) = -1.0 + 2.0 * uGBu - (mux * mux + varx - 1.0);
        }

        const double wi =
            (obs_weight.n_elem == (unsigned)nb) ? obs_weight(ii) : 1.0;
        if (std::isfinite(wi) && wi != 0.0) Glocal += wi * Gi;
    }
    });

    for (int t = 0; t < nt; ++t) G += Gthread[t];

    return Rcpp::List::create(
        Rcpp::Named("log_value") = logv,
        Rcpp::Named("grad_Psi") = G,
        Rcpp::Named("adj_c1") = adjC1,
        Rcpp::Named("adj_c2") = adjC2
    );
}

// [[Rcpp::export]]
Rcpp::List hermite_parallel_info_cpp() {
    const unsigned int hw = std::thread::hardware_concurrency();
#ifdef _OPENMP
    const bool omp_enabled = true;
#else
    const bool omp_enabled = false;
#endif
    return Rcpp::List::create(
        Rcpp::Named("backend") = "std::thread",
        Rcpp::Named("std_thread") = true,
        Rcpp::Named("openmp") = omp_enabled,
        Rcpp::Named("hardware_threads") = (hw > 0 ? (int)hw : NA_INTEGER),
        Rcpp::Named("max_threads") = (hw > 0 ? (int)hw : 1)
    );
}
