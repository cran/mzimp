#include <RcppArmadillo.h>
#include <vector>
#include <cmath>
#include <limits>
#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::plugins(openmp)]]

static inline double ncdf_fast2(const double x) {
    return 0.5 * std::erfc(-x * 0.70710678118654752440);
}
static inline double npdf_fast2(const double x) {
    static const double inv_sqrt_2pi = 0.39894228040143267794;
    return inv_sqrt_2pi * std::exp(-0.5 * x * x);
}

// Exact tensor-GH evaluator specialized to Gaussian count factors, L=0, PM only.
// Improvements over the previous fast kernel:
//   1) PN_j is computed ONCE per block and cached for gradient reuse;
//   2) likelihood uses multiplicative accumulation with periodic row rescaling,
//      avoiding one scalar log() per (observation,node,margin);
//   3) block size is automatically capped to keep PN cache memory bounded.
// [[Rcpp::export]]
Rcpp::List gh_pm_L0_eval_cpp_cached_legacy(
    const Rcpp::List& I_y_list,
    const arma::mat& Psi,
    const arma::mat& Znode,
    const arma::vec& Wnode,
    const arma::vec& t_inner,
    const int d0,
    const int block_size = 500,
    const bool need_grad = true,
    const double tiny = 1e-300,
    const int n_threads = 1,
    const double cache_mb = 256.0,
    const int scale_every = 8
) {
    const int k = Psi.n_rows;
    const int J = Psi.n_cols;
    const int Q = Znode.n_rows;
    if ((int)I_y_list.size() != k) Rcpp::stop("I_y_list length != k");
    if ((int)Znode.n_cols != J) Rcpp::stop("Znode dimension mismatch");
    if ((int)Wnode.n_elem != Q) Rcpp::stop("Wnode length mismatch");
    if ((int)t_inner.n_elem != d0 - 1) Rcpp::stop("t_inner length mismatch");

    std::vector<arma::mat> I(k);
    int n = -1;
    for (int j=0;j<k;++j) {
        I[j] = Rcpp::as<arma::mat>(I_y_list[j]);
        if ((int)I[j].n_cols != d0) Rcpp::stop("I_y matrix has wrong d0");
        if (n < 0) n = I[j].n_rows;
        if ((int)I[j].n_rows != n) Rcpp::stop("I_y row mismatch");
    }

    const int nt = std::max(1, n_threads);
    const double tiny_use = std::max(tiny, std::numeric_limits<double>::min());

    // Build Omega^T and dOmega^T once for current Psi.
    std::vector<arma::mat> OmegaT(k);
    std::vector<arma::mat> dOmegaT;
    if (need_grad) dOmegaT.resize(k);

#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(nt) if(nt > 1)
#endif
    for (int j=0;j<k;++j) {
        const double norm2 = arma::dot(Psi.row(j), Psi.row(j));
        if (!std::isfinite(norm2) || norm2 >= 1.0) continue;
        const double s = std::sqrt(std::max(1e-16, 1.0 - norm2));
        const double s3 = s*s*s;
        arma::mat OmT(d0,Q,arma::fill::zeros);
        arma::mat dOmT;
        if (need_grad) dOmT.zeros(d0,Q*J);
        for (int q=0;q<Q;++q) {
            double eta=0.0;
            for (int r=0;r<J;++r) eta += Znode(q,r)*Psi(j,r);
            std::vector<double> cdf(d0+1,0.0), pdf(d0-1), numer(d0-1);
            cdf[d0]=1.0;
            for (int m=0;m<d0-1;++m) {
                numer[m]=t_inner(m)-eta;
                const double a=numer[m]/s;
                pdf[m]=npdf_fast2(a);
                cdf[m+1]=ncdf_fast2(a);
            }
            for (int m=0;m<d0;++m) OmT(m,q)=std::max(0.0,cdf[m+1]-cdf[m]);
            if (need_grad) {
                for (int r=0;r<J;++r) {
                    double prev=0.0;
                    for (int m=0;m<d0-1;++m) {
                        const double da=-Znode(q,r)/s + numer[m]*Psi(j,r)/s3;
                        const double cur=pdf[m]*da;
                        dOmT(m,q+r*Q)=cur-prev;
                        prev=cur;
                    }
                    dOmT(d0-1,q+r*Q)=-prev;
                }
            }
        }
        OmegaT[j]=std::move(OmT);
        if (need_grad) dOmegaT[j]=std::move(dOmT);
    }

    arma::mat gradPsi(k,J,arma::fill::zeros);
    arma::vec logLi(n,arma::fill::zeros);
    double nll=0.0;

    // Automatically cap block size based on PN cache memory.
    const double bytes_budget = std::max(32.0, cache_mb) * 1024.0 * 1024.0;
    int bs_mem = (int)std::floor(bytes_budget / (8.0 * std::max(1,k) * std::max(1,Q)));
    bs_mem = std::max(8, bs_mem);
    const int bs = std::max(1, std::min(block_size, bs_mem));
    const int se = std::max(1, scale_every);

    for (int i0=0;i0<n;i0+=bs) {
        const int i1=std::min(n-1,i0+bs-1);
        const int nb=i1-i0+1;
        std::vector<arma::mat> PN(k);

        // Compute PN once per margin and cache it.
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(nt) if(nt > 1)
#endif
        for (int j=0;j<k;++j) {
            arma::mat Ib = I[j].rows(i0,i1);
            PN[j] = ((double)d0) * (Ib * OmegaT[j]);
            PN[j].transform([&](double v){ return (v < tiny_use || !std::isfinite(v)) ? tiny_use : v; });
        }

        // Multiplicative integrand with periodic row scaling.
        arma::mat T(nb,Q,arma::fill::ones);
        arma::vec logscale(nb,arma::fill::zeros);
        for (int j=0;j<k;++j) {
            T %= PN[j];
            if (((j+1)%se)==0 || j==k-1) {
                arma::vec mx = arma::max(T,1);
                for (int ii=0;ii<nb;++ii) {
                    double m = mx(ii);
                    if (!std::isfinite(m) || m <= 0.0) m = tiny_use;
                    T.row(ii) /= m;
                    logscale(ii) += std::log(m);
                }
            }
        }

        // Node posterior weights alpha and observation log likelihoods.
        arma::mat alpha = T;
        for (int q=0;q<Q;++q) alpha.col(q) *= Wnode(q);
        arma::vec denom = arma::sum(alpha,1);
        for (int ii=0;ii<nb;++ii) {
            double ss=denom(ii);
            if (!std::isfinite(ss) || ss<=0.0) Rcpp::stop("GH normalization failed");
            const double lli=logscale(ii)+std::log(ss);
            logLi(i0+ii)=lli;
            nll-=lli;
            alpha.row(ii)/=ss;
        }

        if (need_grad) {
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(nt) if(nt > 1)
#endif
            for (int j=0;j<k;++j) {
                arma::mat Ib=I[j].rows(i0,i1);
                arma::mat dPNall=((double)d0)*(Ib*dOmegaT[j]);
                for (int r=0;r<J;++r) {
                    const int c0=r*Q;
                    double g=0.0;
                    for (int ii=0;ii<nb;++ii) {
                        // Armadillo rows are strided in column-major storage; use indexing.
                        for (int q=0;q<Q;++q)
                            g -= alpha(ii,q) * dPNall(ii,c0+q) / PN[j](ii,q);
                    }
                    gradPsi(j,r)+=g;
                }
            }
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
        Rcpp::Named("Q")=Q,
        Rcpp::Named("threads")=nt,
        Rcpp::Named("openmp")=omp_enabled,
        Rcpp::Named("effective_block_size")=bs
    );
}


// Observation-block OpenMP evaluator. Used only when n_threads > 1.
// Each OpenMP thread owns complete observation blocks and performs PN,
// likelihood, posterior weights, and gradient for those blocks.
// [[Rcpp::export]]
Rcpp::List gh_pm_L0_eval_cpp_cached_obsblock(
    const Rcpp::List& I_y_list,
    const arma::mat& Psi,
    const arma::mat& Znode,
    const arma::vec& Wnode,
    const arma::vec& t_inner,
    const int d0,
    const int block_size = 500,
    const bool need_grad = true,
    const double tiny = 1e-300,
    const int n_threads = 2,
    const double cache_mb = 256.0,
    const int scale_every = 8
) {
    const int k = Psi.n_rows;
    const int J = Psi.n_cols;
    const int Q = Znode.n_rows;
    if ((int)I_y_list.size() != k) Rcpp::stop("I_y_list length != k");
    if ((int)Znode.n_cols != J) Rcpp::stop("Znode dimension mismatch");
    if ((int)Wnode.n_elem != Q) Rcpp::stop("Wnode length mismatch");
    if ((int)t_inner.n_elem != d0 - 1) Rcpp::stop("t_inner length mismatch");

    std::vector<arma::mat> I(k);
    int n = -1;
    for (int j=0; j<k; ++j) {
        I[j] = Rcpp::as<arma::mat>(I_y_list[j]);
        if ((int)I[j].n_cols != d0) Rcpp::stop("I_y matrix has wrong d0");
        if (n < 0) n = I[j].n_rows;
        if ((int)I[j].n_rows != n) Rcpp::stop("I_y row mismatch");
    }

    const int nt = std::max(1, n_threads);
    const double tiny_use = std::max(tiny, std::numeric_limits<double>::min());
    const int se = std::max(1, scale_every);

    // Omega depends on Psi but not observations: build once per evaluation.
    // Keep this serial in the new path to avoid nested/competing parallel regions.
    std::vector<arma::mat> OmegaT(k);
    std::vector<arma::mat> dOmegaT;
    if (need_grad) dOmegaT.resize(k);
    for (int j=0; j<k; ++j) {
        const double norm2 = arma::dot(Psi.row(j), Psi.row(j));
        if (!std::isfinite(norm2) || norm2 >= 1.0) continue;
        const double ss = std::sqrt(std::max(1e-16, 1.0 - norm2));
        const double s3 = ss*ss*ss;
        arma::mat OmT(d0,Q,arma::fill::zeros);
        arma::mat dOmT;
        if (need_grad) dOmT.zeros(d0,Q*J);
        for (int q=0; q<Q; ++q) {
            double eta=0.0;
            for (int r=0; r<J; ++r) eta += Znode(q,r)*Psi(j,r);
            std::vector<double> cdf(d0+1,0.0), pdf(d0-1), numer(d0-1);
            cdf[d0]=1.0;
            for (int m=0; m<d0-1; ++m) {
                numer[m]=t_inner(m)-eta;
                const double a=numer[m]/ss;
                pdf[m]=npdf_fast2(a);
                cdf[m+1]=ncdf_fast2(a);
            }
            for (int m=0; m<d0; ++m) OmT(m,q)=std::max(0.0,cdf[m+1]-cdf[m]);
            if (need_grad) {
                for (int r=0; r<J; ++r) {
                    double prev=0.0;
                    for (int m=0; m<d0-1; ++m) {
                        const double da=-Znode(q,r)/ss + numer[m]*Psi(j,r)/s3;
                        const double cur=pdf[m]*da;
                        dOmT(m,q+r*Q)=cur-prev;
                        prev=cur;
                    }
                    dOmT(d0-1,q+r*Q)=-prev;
                }
            }
        }
        OmegaT[j]=std::move(OmT);
        if (need_grad) dOmegaT[j]=std::move(dOmT);
    }

    arma::mat gradPsi(k,J,arma::fill::zeros);
    arma::vec logLi(n,arma::fill::zeros);
    double nll=0.0;

    // cache_mb is a TOTAL PN-cache budget. Divide it across active threads.
    const double total_bytes = std::max(32.0, cache_mb) * 1024.0 * 1024.0;
    const double bytes_per_thread = total_bytes / (double)nt;
    int bs_mem = (int)std::floor(bytes_per_thread /
        (8.0 * std::max(1,k) * std::max(1,Q)));
    bs_mem = std::max(1, bs_mem);

    // For parallel mode, target roughly one primary block per thread when possible.
    // block_size remains an upper bound supplied by R/user.
    const int bs_balanced = std::max(1, (n + nt - 1) / nt);
    const int bs = std::max(1, std::min(block_size, std::min(bs_mem, bs_balanced)));
    const int nblocks = (n + bs - 1) / bs;

    std::vector<double> nll_thread(nt, 0.0);
    std::vector<arma::mat> grad_thread;
    if (need_grad) {
        grad_thread.reserve(nt);
        for (int t=0; t<nt; ++t) grad_thread.emplace_back(k,J,arma::fill::zeros);
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
        arma::mat grad_local;
        if (need_grad) grad_local.zeros(k,J);

#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
        for (int b=0; b<nblocks; ++b) {
            const int i0 = b*bs;
            const int i1 = std::min(n-1, i0+bs-1);
            const int nb = i1-i0+1;
            std::vector<arma::mat> PN(k);

            // All margins for this observation block are owned by this thread.
            for (int j=0; j<k; ++j) {
                arma::mat Ib = I[j].rows(i0,i1);
                PN[j] = ((double)d0) * (Ib * OmegaT[j]);
                PN[j].transform([&](double v){
                    return (v < tiny_use || !std::isfinite(v)) ? tiny_use : v;
                });
            }

            arma::mat T(nb,Q,arma::fill::ones);
            arma::vec logscale(nb,arma::fill::zeros);
            for (int j=0; j<k; ++j) {
                T %= PN[j];
                if (((j+1)%se)==0 || j==k-1) {
                    arma::vec mx = arma::max(T,1);
                    for (int ii=0; ii<nb; ++ii) {
                        double mm = mx(ii);
                        if (!std::isfinite(mm) || mm <= 0.0) mm = tiny_use;
                        T.row(ii) /= mm;
                        logscale(ii) += std::log(mm);
                    }
                }
            }

            arma::mat alpha = T;
            for (int q=0; q<Q; ++q) alpha.col(q) *= Wnode(q);
            arma::vec denom = arma::sum(alpha,1);
            bool block_ok = true;
            for (int ii=0; ii<nb; ++ii) {
                const double den = denom(ii);
                if (!std::isfinite(den) || den <= 0.0) {
                    block_ok = false;
                    continue;
                }
                const double lli = logscale(ii) + std::log(den);
                logLi(i0+ii) = lli; // disjoint rows across blocks: thread-safe
                nll_local -= lli;
                alpha.row(ii) /= den;
            }
            if (!block_ok) {
#ifdef _OPENMP
#pragma omp atomic write
#endif
                failed = 1;
                continue;
            }

            if (need_grad) {
                for (int j=0; j<k; ++j) {
                    arma::mat Ib = I[j].rows(i0,i1);
                    arma::mat dPNall = ((double)d0) * (Ib * dOmegaT[j]);
                    for (int r=0; r<J; ++r) {
                        const int c0 = r*Q;
                        double g=0.0;
                        for (int ii=0; ii<nb; ++ii) {
                            for (int q=0; q<Q; ++q)
                                g -= alpha(ii,q) * dPNall(ii,c0+q) / PN[j](ii,q);
                        }
                        grad_local(j,r) += g;
                    }
                }
            }
        }

        nll_thread[tid] = nll_local;
        if (need_grad) grad_thread[tid] = std::move(grad_local);
    }

    if (failed) Rcpp::stop("GH normalization failed");
    for (int t=0; t<nt; ++t) {
        nll += nll_thread[t];
        if (need_grad) gradPsi += grad_thread[t];
    }

    bool omp_enabled=false;
#ifdef _OPENMP
    omp_enabled=true;
#endif
    return Rcpp::List::create(
        Rcpp::Named("nll")=nll,
        Rcpp::Named("logLi")=logLi,
        Rcpp::Named("grad_Psi")=gradPsi,
        Rcpp::Named("Q")=Q,
        Rcpp::Named("threads")=nt,
        Rcpp::Named("openmp")=omp_enabled,
        Rcpp::Named("effective_block_size")=bs,
        Rcpp::Named("n_blocks")=nblocks,
        Rcpp::Named("parallel_axis")="observation_blocks"
    );
}

// Public compatibility wrapper:
// n_threads == 1 -> EXACT legacy path, unchanged.
// n_threads > 1  -> new observation-block OpenMP path.
// [[Rcpp::export]]
Rcpp::List gh_pm_L0_eval_cpp_cached(
    const Rcpp::List& I_y_list,
    const arma::mat& Psi,
    const arma::mat& Znode,
    const arma::vec& Wnode,
    const arma::vec& t_inner,
    const int d0,
    const int block_size = 500,
    const bool need_grad = true,
    const double tiny = 1e-300,
    const int n_threads = 1,
    const double cache_mb = 256.0,
    const int scale_every = 8
) {
    if (n_threads <= 1) {
        return gh_pm_L0_eval_cpp_cached_legacy(
            I_y_list, Psi, Znode, Wnode, t_inner, d0,
            block_size, need_grad, tiny, 1, cache_mb, scale_every
        );
    }
    return gh_pm_L0_eval_cpp_cached_obsblock(
        I_y_list, Psi, Znode, Wnode, t_inner, d0,
        block_size, need_grad, tiny, n_threads, cache_mb, scale_every
    );
}

// [[Rcpp::export]]
Rcpp::List gh_fast_cached_parallel_info_cpp() {
#ifdef _OPENMP
    return Rcpp::List::create(Rcpp::Named("openmp")=true,Rcpp::Named("max_threads")=omp_get_max_threads());
#else
    return Rcpp::List::create(Rcpp::Named("openmp")=false,Rcpp::Named("max_threads")=1);
#endif
}
