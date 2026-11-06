#include <RcppArmadillo.h>
#include <vector>
#include <cmath>
#include <limits>
#include <algorithm>
#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::plugins(openmp)]]

static inline double ghL1_cdf(const double x) {
    return 0.5 * std::erfc(-x * 0.70710678118654752440);
}
static inline double ghL1_pdf(const double x) {
    static const double inv_sqrt_2pi = 0.39894228040143267794;
    return inv_sqrt_2pi * std::exp(-0.5 * x * x);
}
static inline double ghL1_pos(const double x, const double tiny) {
    return (!std::isfinite(x) || x < tiny) ? tiny : x;
}

// Exact tensor-GH evaluator for Gaussian count factors + Gaussian L=1 zero factor.
//
// IMPORTANT DESIGN:
//   * This backend is intended for the active L=1 structural-zero factor case.
//   * Parallelism is over OBSERVATION BLOCKS, exactly like the new L=0 backend.
//   * Each OpenMP thread owns complete observation blocks and computes:
//       PN_y / PN_0 -> PB -> q x b contraction -> logLik -> gradients.
//   * No R-level per-block parallelism and no nb x (QN*QB) R object is materialized.
//   * qpb_list contains qnorm(pb) for ZIPM margins only, in ZIPM order.
//   * zi_pos is length k, with 0 for PM and 1..nzi for ZIPM margins.
//
// Returned gradients are with respect to the constrained loading matrices Psi and
// Rho, matching the old gh_L1_fused_cpp() convention.  R can apply the existing
// chain_grad_loading_to_raw() afterward.
//
// [[Rcpp::export]]
Rcpp::List gh_L1_eval_cpp_cached_obsblock(
    const Rcpp::List& I_y_list,
    const Rcpp::List& I_0_list,
    const Rcpp::List& qpb_list,
    const Rcpp::List& y0_list,
    const Rcpp::IntegerVector& zi_pos,
    const arma::mat& Psi,
    const arma::vec& Rho,
    const arma::mat& Znode_N,
    const arma::vec& Wnode_N,
    const arma::vec& Znode_B,
    const arma::vec& Wnode_B,
    const arma::vec& t_inner,
    const int d0,
    const int block_size = 500,
    const bool need_grad = true,
    const double tiny = 1e-300,
    const int n_threads = 2,
    const double cache_mb = 256.0
) {
    const int k  = Psi.n_rows;
    const int J  = Psi.n_cols;
    const int QN = Znode_N.n_rows;
    const int QB = Znode_B.n_elem;
    const int nzi = Rho.n_elem;

    if ((int)I_y_list.size() != k) Rcpp::stop("I_y_list length != k");
    if ((int)I_0_list.size() != k) Rcpp::stop("I_0_list length != k");
    if ((int)y0_list.size() != k) Rcpp::stop("y0_list length != k");
    if ((int)zi_pos.size() != k) Rcpp::stop("zi_pos length != k");
    if ((int)qpb_list.size() != nzi) Rcpp::stop("qpb_list length != length(Rho)");
    if ((int)Znode_N.n_cols != J) Rcpp::stop("Znode_N dimension mismatch");
    if ((int)Wnode_N.n_elem != QN) Rcpp::stop("Wnode_N length mismatch");
    if ((int)Wnode_B.n_elem != QB) Rcpp::stop("Wnode_B length mismatch");
    if ((int)t_inner.n_elem != d0 - 1) Rcpp::stop("t_inner length mismatch");
    if (d0 < 2) Rcpp::stop("d0 must be >= 2");
    if (QN < 1 || QB < 1) Rcpp::stop("GH node grids must be non-empty");

    // Copy R inputs once before entering OpenMP.  No R API is touched in workers.
    std::vector<arma::mat> Iy(k), I0(k);
    std::vector<arma::uvec> y0(k);
    std::vector<arma::vec> qpb(nzi);
    int n = -1;

    for (int j=0; j<k; ++j) {
        Iy[j] = Rcpp::as<arma::mat>(I_y_list[j]);
        if ((int)Iy[j].n_cols != d0) Rcpp::stop("I_y matrix has wrong d0");
        if (n < 0) n = Iy[j].n_rows;
        if ((int)Iy[j].n_rows != n) Rcpp::stop("I_y row mismatch");

        Rcpp::LogicalVector yj = y0_list[j];
        if ((int)yj.size() != n) Rcpp::stop("y0 length mismatch");
        y0[j].set_size(n);
        for (int i=0; i<n; ++i) y0[j](i) = (yj[i] == TRUE) ? 1u : 0u;

        const int s1 = zi_pos[j];
        if (s1 < 0 || s1 > nzi) Rcpp::stop("Invalid zi_pos entry");
        if (s1 > 0) {
            I0[j] = Rcpp::as<arma::mat>(I_0_list[j]);
            if ((int)I0[j].n_rows != n || (int)I0[j].n_cols != d0)
                Rcpp::stop("I_0 matrix dimension mismatch");
        }
    }
    for (int s=0; s<nzi; ++s) {
        qpb[s] = Rcpp::as<arma::vec>(qpb_list[s]);
        if ((int)qpb[s].n_elem != n) Rcpp::stop("qpb vector length mismatch");
    }

    const int nt = std::max(1, n_threads);
    const double tiny_use = std::max(tiny, std::numeric_limits<double>::min());

    // ---------------------------------------------------------------------
    // Count-factor cell probabilities Omega_N and derivatives.
    // Depends on parameters/nodes, not observations: build once.
    // ---------------------------------------------------------------------
    std::vector<arma::mat> OmegaT(k);
    std::vector<arma::mat> dOmegaT;
    if (need_grad) dOmegaT.resize(k);

    for (int j=0; j<k; ++j) {
        const double norm2 = arma::dot(Psi.row(j), Psi.row(j));
        if (!std::isfinite(norm2) || norm2 >= 1.0)
            Rcpp::stop("Invalid count-factor loading matrix");
        const double sig = std::sqrt(std::max(1e-16, 1.0 - norm2));
        const double sig3 = sig*sig*sig;

        arma::mat OmT(d0, QN, arma::fill::zeros);
        arma::mat dOmT;
        if (need_grad) dOmT.zeros(d0, QN*J);

        for (int q=0; q<QN; ++q) {
            double eta = 0.0;
            for (int r=0; r<J; ++r) eta += Znode_N(q,r) * Psi(j,r);

            std::vector<double> cdf(d0+1, 0.0), pdf(d0-1), numer(d0-1);
            cdf[d0] = 1.0;
            for (int m=0; m<d0-1; ++m) {
                numer[m] = t_inner(m) - eta;
                const double a = numer[m] / sig;
                pdf[m] = ghL1_pdf(a);
                cdf[m+1] = ghL1_cdf(a);
            }
            for (int m=0; m<d0; ++m)
                OmT(m,q) = std::max(0.0, cdf[m+1] - cdf[m]);

            if (need_grad) {
                for (int r=0; r<J; ++r) {
                    double prev = 0.0;
                    for (int m=0; m<d0-1; ++m) {
                        const double da = -Znode_N(q,r)/sig + numer[m]*Psi(j,r)/sig3;
                        const double cur = pdf[m] * da;
                        dOmT(m, q + r*QN) = cur - prev;
                        prev = cur;
                    }
                    dOmT(d0-1, q + r*QN) = -prev;
                }
            }
        }
        OmegaT[j] = std::move(OmT);
        if (need_grad) dOmegaT[j] = std::move(dOmT);
    }

    // ---------------------------------------------------------------------
    // Zero-factor constants for L=1 Gaussian factor.
    // PB_i,s,b = Phi((qpb_i,s - rho_s*z_b)/sqrt(1-rho_s^2)).
    // ---------------------------------------------------------------------
    arma::vec sigB(nzi), sigB3(nzi);
    for (int s=0; s<nzi; ++s) {
        const double r = Rho(s);
        if (!std::isfinite(r) || r*r >= 1.0)
            Rcpp::stop("Invalid zero-factor loading vector");
        sigB(s) = std::sqrt(std::max(1e-16, 1.0-r*r));
        sigB3(s) = sigB(s)*sigB(s)*sigB(s);
    }

    arma::mat gradPsi(k, J, arma::fill::zeros);
    arma::vec gradRho(nzi, arma::fill::zeros);
    arma::vec logLi(n, arma::fill::zeros);
    double nll = 0.0;

    // Memory estimate per observation includes PN_y for all margins and PN_0
    // for ZIPM margins, plus PB and posterior workspace.  cache_mb is total.
    const double total_bytes = std::max(32.0, cache_mb) * 1024.0 * 1024.0;
    const double bytes_per_thread = total_bytes / (double)nt;
    const double per_obs_doubles =
        (double)k * QN + (double)nzi * QN + (double)nzi * QB +
        2.0*(double)QN*(double)QB + 2.0*(double)(QN+QB) + 16.0;
    int bs_mem = (int)std::floor(bytes_per_thread /
        (8.0 * std::max(1.0, per_obs_doubles)));
    bs_mem = std::max(1, bs_mem);
    const int bs_balanced = std::max(1, (n + nt - 1) / nt);
    const int bs = std::max(1, std::min(block_size, std::min(bs_mem, bs_balanced)));
    const int nblocks = (n + bs - 1) / bs;

    std::vector<double> nll_thread(nt, 0.0);
    std::vector<arma::mat> gPsi_thread;
    std::vector<arma::vec> gRho_thread;
    if (need_grad) {
        gPsi_thread.reserve(nt);
        gRho_thread.reserve(nt);
        for (int t=0; t<nt; ++t) {
            gPsi_thread.emplace_back(k, J, arma::fill::zeros);
            gRho_thread.emplace_back(nzi, arma::fill::zeros);
        }
    }
    int failed = 0;

#ifdef _OPENMP
#pragma omp parallel num_threads(nt) if(nt > 1)
#endif
    {
        int tid = 0;
#ifdef _OPENMP
        tid = omp_get_thread_num();
#endif
        double nll_local = 0.0;
        arma::mat gPsi_local;
        arma::vec gRho_local;
        if (need_grad) {
            gPsi_local.zeros(k,J);
            gRho_local.zeros(nzi);
        }

#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
        for (int bb=0; bb<nblocks; ++bb) {
            const int i0 = bb*bs;
            const int i1 = std::min(n-1, i0+bs-1);
            const int nb = i1-i0+1;

            // PN caches: all margins for this observation block belong to this thread.
            std::vector<arma::mat> PNy(k), PN0(k);
            for (int j=0; j<k; ++j) {
                arma::mat Iyb = Iy[j].rows(i0,i1);
                PNy[j] = ((double)d0) * (Iyb * OmegaT[j]);
                PNy[j].transform([&](double v){ return ghL1_pos(v, tiny_use); });
                if (zi_pos[j] > 0) {
                    arma::mat I0b = I0[j].rows(i0,i1);
                    PN0[j] = ((double)d0) * (I0b * OmegaT[j]);
                    PN0[j].transform([&](double v){ return ghL1_pos(v, tiny_use); });
                }
            }

            // PB caches for ZIPM margins: nb x QB.
            std::vector<arma::mat> PB(nzi), dPB;
            if (need_grad) dPB.resize(nzi);
            for (int s=0; s<nzi; ++s) {
                PB[s].set_size(nb,QB);
                if (need_grad) dPB[s].set_size(nb,QB);
                const double rho = Rho(s);
                const double sig = sigB(s);
                const double sig3 = sigB3(s);
                for (int b=0; b<QB; ++b) {
                    const double z = Znode_B(b);
                    for (int ii=0; ii<nb; ++ii) {
                        const double A = qpb[s](i0+ii) - rho*z;
                        const double av = A/sig;
                        PB[s](ii,b) = ghL1_cdf(av);
                        if (need_grad) {
                            const double da = A*rho/sig3 - z/sig;
                            dPB[s](ii,b) = ghL1_pdf(av) * da;
                        }
                    }
                }
            }

            // Per-observation exact q x b contraction.
            std::vector<double> Aq(QN), Bb(QB);
            std::vector<double> logK((size_t)QN*QB);
            std::vector<double> alpha((size_t)QN*QB);
            std::vector<double> alphaN(QN), alphaB(QB);
            std::vector<int> zero_js;
            zero_js.reserve(k);

            for (int ii=0; ii<nb; ++ii) {
                const int ig = i0+ii;
                for (int q=0; q<QN; ++q) Aq[q] = std::log(ghL1_pos(Wnode_N(q), tiny_use));
                for (int b=0; b<QB; ++b) Bb[b] = std::log(ghL1_pos(Wnode_B(b), tiny_use));
                zero_js.clear();

                for (int j=0; j<k; ++j) {
                    const int s1 = zi_pos[j];
                    if (s1 == 0) {
                        for (int q=0; q<QN; ++q)
                            Aq[q] += std::log(PNy[j](ii,q));
                    } else if (y0[j](ig) != 0u) {
                        zero_js.push_back(j);
                    } else {
                        const int s = s1-1;
                        for (int q=0; q<QN; ++q)
                            Aq[q] += std::log(PNy[j](ii,q));
                        for (int b=0; b<QB; ++b)
                            Bb[b] += std::log(ghL1_pos(PB[s](ii,b), tiny_use));
                    }
                }

                double mx = -std::numeric_limits<double>::infinity();
                for (int q=0; q<QN; ++q) {
                    for (int b=0; b<QB; ++b) {
                        double lk = Aq[q] + Bb[b];
                        for (size_t zz=0; zz<zero_js.size(); ++zz) {
                            const int j = zero_js[zz];
                            const int s = zi_pos[j]-1;
                            const double pb = PB[s](ii,b);
                            const double p0 = PN0[j](ii,q);
                            const double T = ghL1_pos(1.0 + pb*(p0-1.0), tiny_use);
                            lk += std::log(T);
                        }
                        const size_t id = (size_t)q*QB+b;
                        logK[id] = lk;
                        if (lk > mx) mx = lk;
                    }
                }

                if (!std::isfinite(mx)) {
#ifdef _OPENMP
#pragma omp atomic write
#endif
                    failed = 1;
                    continue;
                }
                double den = 0.0;
                for (size_t tt=0; tt<logK.size(); ++tt) den += std::exp(logK[tt]-mx);
                if (!(den > 0.0) || !std::isfinite(den)) {
#ifdef _OPENMP
#pragma omp atomic write
#endif
                    failed = 1;
                    continue;
                }
                const double lli = mx + std::log(den);
                logLi(ig) = lli; // disjoint observation rows
                nll_local -= lli;

                if (!need_grad) continue;

                std::fill(alphaN.begin(), alphaN.end(), 0.0);
                std::fill(alphaB.begin(), alphaB.end(), 0.0);
                for (int q=0; q<QN; ++q) {
                    for (int b=0; b<QB; ++b) {
                        const size_t id=(size_t)q*QB+b;
                        const double a=std::exp(logK[id]-lli);
                        alpha[id]=a;
                        alphaN[q]+=a;
                        alphaB[b]+=a;
                    }
                }

                for (int j=0; j<k; ++j) {
                    const int s1=zi_pos[j];
                    arma::rowvec Iyrow = Iy[j].row(ig);
                    arma::rowvec I0row;
                    if (s1>0) I0row = I0[j].row(ig);

                    if (s1==0 || y0[j](ig)==0u) {
                        // PM or positive ZIPM: count score depends only on q.
                        if (need_grad) {
                            arma::rowvec dPNall = ((double)d0) * (Iyrow * dOmegaT[j]);
                            for (int r=0; r<J; ++r) {
                                double sc=0.0;
                                const int off=r*QN;
                                for (int q=0; q<QN; ++q)
                                    sc += alphaN[q] * dPNall(off+q) / PNy[j](ii,q);
                                gPsi_local(j,r) -= sc;
                            }
                        }
                        if (s1>0) {
                            const int s=s1-1;
                            double sr=0.0;
                            for (int b=0; b<QB; ++b) {
                                const double pb=ghL1_pos(PB[s](ii,b), tiny_use);
                                sr += alphaB[b] * dPB[s](ii,b) / pb;
                            }
                            gRho_local(s) -= sr;
                        }
                    } else {
                        // Observed ZIPM zero: genuinely coupled q-b score.
                        const int s=s1-1;
                        arma::rowvec dPN0all = ((double)d0) * (I0row * dOmegaT[j]);
                        for (int r=0; r<J; ++r) {
                            double sc=0.0;
                            const int off=r*QN;
                            for (int q=0; q<QN; ++q) {
                                const double p0=PN0[j](ii,q);
                                const double dp0=dPN0all(off+q);
                                for (int b=0; b<QB; ++b) {
                                    const double pb=PB[s](ii,b);
                                    const double T=ghL1_pos(1.0+pb*(p0-1.0),tiny_use);
                                    sc += alpha[(size_t)q*QB+b] * (pb*dp0/T);
                                }
                            }
                            gPsi_local(j,r) -= sc;
                        }
                        double sr=0.0;
                        for (int q=0; q<QN; ++q) {
                            const double p0=PN0[j](ii,q);
                            for (int b=0; b<QB; ++b) {
                                const double pb=PB[s](ii,b);
                                const double T=ghL1_pos(1.0+pb*(p0-1.0),tiny_use);
                                sr += alpha[(size_t)q*QB+b] * dPB[s](ii,b)*(p0-1.0)/T;
                            }
                        }
                        gRho_local(s) -= sr;
                    }
                }
            }
        }

        nll_thread[tid] = nll_local;
        if (need_grad) {
            gPsi_thread[tid] = std::move(gPsi_local);
            gRho_thread[tid] = std::move(gRho_local);
        }
    }

    if (failed) Rcpp::stop("L=1 GH normalization failed");
    for (int t=0; t<nt; ++t) {
        nll += nll_thread[t];
        if (need_grad) {
            gradPsi += gPsi_thread[t];
            gradRho += gRho_thread[t];
        }
    }

    bool omp_enabled=false;
#ifdef _OPENMP
    omp_enabled=true;
#endif
    return Rcpp::List::create(
        Rcpp::Named("nll")=nll,
        Rcpp::Named("logLi")=logLi,
        Rcpp::Named("grad_Psi")=gradPsi,
        Rcpp::Named("grad_Rho")=gradRho,
        Rcpp::Named("QN")=QN,
        Rcpp::Named("QB")=QB,
        Rcpp::Named("threads")=nt,
        Rcpp::Named("openmp")=omp_enabled,
        Rcpp::Named("effective_block_size")=bs,
        Rcpp::Named("n_blocks")=nblocks,
        Rcpp::Named("parallel_axis")="observation_blocks"
    );
}

// Convenience wrapper.  The R layer can choose to call this only when
// parallel > 1 and keep the existing legacy R + gh_L1_fused_cpp path for
// parallel == 1.  If called with n_threads=1, this same backend runs serially.
// [[Rcpp::export]]
Rcpp::List gh_L1_eval_cpp_cached(
    const Rcpp::List& I_y_list,
    const Rcpp::List& I_0_list,
    const Rcpp::List& qpb_list,
    const Rcpp::List& y0_list,
    const Rcpp::IntegerVector& zi_pos,
    const arma::mat& Psi,
    const arma::vec& Rho,
    const arma::mat& Znode_N,
    const arma::vec& Wnode_N,
    const arma::vec& Znode_B,
    const arma::vec& Wnode_B,
    const arma::vec& t_inner,
    const int d0,
    const int block_size = 500,
    const bool need_grad = true,
    const double tiny = 1e-300,
    const int n_threads = 1,
    const double cache_mb = 256.0
) {
    return gh_L1_eval_cpp_cached_obsblock(
        I_y_list, I_0_list, qpb_list, y0_list, zi_pos,
        Psi, Rho, Znode_N, Wnode_N, Znode_B, Wnode_B,
        t_inner, d0, block_size, need_grad, tiny, n_threads, cache_mb
    );
}

// [[Rcpp::export]]
Rcpp::List gh_L1_parallel_info_cpp() {
#ifdef _OPENMP
    return Rcpp::List::create(
        Rcpp::Named("openmp")=true,
        Rcpp::Named("max_threads")=omp_get_max_threads(),
        Rcpp::Named("parallel_axis")="observation_blocks"
    );
#else
    return Rcpp::List::create(
        Rcpp::Named("openmp")=false,
        Rcpp::Named("max_threads")=1,
        Rcpp::Named("parallel_axis")="serial_fallback"
    );
#endif
}
