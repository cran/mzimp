#include <Rcpp.h>
#include <cmath>
#include <thread>
#include <vector>
#include <algorithm>
#include <limits>
// [[Rcpp::plugins(cpp11)]]
using namespace Rcpp;

namespace {
inline double Phi(double x){ return 0.5*std::erfc(-x*M_SQRT1_2); }
inline double clampu(double u){ return std::min(1.0-1e-15,std::max(1e-15,u)); }
inline double logadd(double a,double b){ if(!std::isfinite(a)) return b; if(!std::isfinite(b)) return a; double m=std::max(a,b); return m+std::log(std::exp(a-m)+std::exp(b-m)); }
inline int nth(int req,int n){ return std::max(1,std::min(std::max(1,req),n)); }
inline double qnorm_safe(double p,double mu,double sd){ return R::qnorm5(clampu(p),mu,sd,1,0); }
}

// -----------------------------------------------------------------------------
// Sequential conditional (Genz-style) RQMC likelihood for Gaussian factor ZIPM.
// Supports:
//   * L=0: PM + independent ZIPM margins (pb_list supplies active probability),
//   * L=1: only ZIPM margins indexed by zi_position use the zero factor Rho.
// Non-ZI margins are handled as ordinary count margins.
//
// Convention inherited from the Hermite code:
//   PB = P(count component is active).
// Thus for y=0: P = (1-PB) + PB * P_count(Y*=0).
//
// U_count and U_zero are Q x k randomized low-discrepancy uniforms.  U_zero is
// ignored for non-ZI margins and when L=0 except for selecting the zero mixture
// branch.  Parallelization is over observations; each thread owns all particle
// state, so no R API is called from worker threads except R::qnorm5 (thread-safe
// math routine in Rmath in normal package builds).
// -----------------------------------------------------------------------------
// [[Rcpp::export]]
Rcpp::NumericVector genz_conditional_zipm_qmc_cpp(
    const Rcpp::NumericMatrix& U_count,
    const Rcpp::NumericMatrix& U_zero,
    const Rcpp::NumericMatrix& Psi,
    const Rcpp::NumericMatrix& Rho,
    const Rcpp::IntegerVector& zi_position,
    const Rcpp::NumericVector& t_inner,
    const Rcpp::List& I_y_list,
    const Rcpp::List& I_0_list,
    const Rcpp::LogicalVector& zero_inflation,
    const Rcpp::List& pb_list,
    const Rcpp::List& qpb_list,
    const Rcpp::List& y0_list,
    const double tiny=1e-300,
    const int n_threads=1) {

  const int Q=U_count.nrow(), k=U_count.ncol(), J=Psi.ncol();
  if(U_zero.nrow()!=Q || U_zero.ncol()!=k) stop("U_zero must have same Q x k dimension as U_count.");
  if(Psi.nrow()!=k || I_y_list.size()!=k || I_0_list.size()!=k || zero_inflation.size()!=k ||
     pb_list.size()!=k || qpb_list.size()!=k || y0_list.size()!=k || zi_position.size()!=k) stop("dimension mismatch");
  const int L=Rho.ncol();
  if(L<0 || L>1) stop("This Genz kernel supports only L=0 or L=1.");
  NumericMatrix A0=I_y_list[0]; const int n=A0.nrow(), d=A0.ncol();
  if(t_inner.size()!=std::max(0,d-1)) stop("t_inner mismatch");

  std::vector<double> edge(d + 1);
  edge[0] = -std::numeric_limits<double>::infinity();

  for (int m = 1; m < d; ++m) {
    edge[m] = t_inner[m - 1];
  }

  edge[d] = std::numeric_limits<double>::infinity();  

  // Copy all R inputs needed by workers to plain C++ memory.
  std::vector<double> psi((size_t)k*J); std::copy(Psi.begin(),Psi.end(),psi.begin());
  std::vector<double> rho;
  const int nzi=Rho.nrow(); if(L==1){ rho.resize(nzi); for(int s=0;s<nzi;++s) rho[s]=Rho(s,0); }
  std::vector<int> zip(k), zipos(k), yzero((size_t)n*k,0);
  std::vector< std::vector<double> > Iy(k), I0(k), pb(k), qpb(k);
  for(int j=0;j<k;++j){
    NumericMatrix Ay=I_y_list[j]; if(Ay.nrow()!=n || Ay.ncol()!=d) stop("I_y dimensions differ");
    Iy[j].assign(Ay.begin(),Ay.end()); zip[j]=(zero_inflation[j]==TRUE); zipos[j]=zi_position[j]-1;
    if(zip[j]){
      NumericMatrix Az=I_0_list[j]; LogicalVector yz=y0_list[j];
      if(Az.nrow()!=n || Az.ncol()!=d || yz.size()!=n) stop("ZIPM dimensions differ");
      I0[j].assign(Az.begin(),Az.end()); for(int i=0;i<n;++i) yzero[(size_t)i*k+j]=(yz[i]==TRUE);
      if(L==0){ NumericVector pp=pb_list[j]; if(pp.size()!=n) stop("pb_list length mismatch"); pb[j].assign(pp.begin(),pp.end()); }
      else { NumericVector qq=qpb_list[j]; if(qq.size()!=n || zipos[j]<0 || zipos[j]>=nzi) stop("qpb/zi_position mismatch"); qpb[j].assign(qq.begin(),qq.end()); }
    }
  }

  // Count residual variances and deterministic Kalman gains/covariance sequence.
  std::vector<double> Dc(k), sc(k), Kc((size_t)k*J,0.0), Vc((size_t)J*J,0.0), vp(J);
  for(int a=0;a<J;++a) Vc[(size_t)a*J+a]=1.0;
  for(int j=0;j<k;++j){ double ss=0; for(int a=0;a<J;++a){ double v=psi[j+(size_t)k*a]; ss+=v*v; } Dc[j]=std::max(1e-12,1.0-ss);
    for(int a=0;a<J;++a){ double z=0; for(int b=0;b<J;++b) z+=Vc[(size_t)a*J+b]*psi[j+(size_t)k*b]; vp[a]=z; }
    double s2=Dc[j]; for(int a=0;a<J;++a) s2+=psi[j+(size_t)k*a]*vp[a]; s2=std::max(1e-12,s2); sc[j]=std::sqrt(s2);
    for(int a=0;a<J;++a) Kc[(size_t)j*J+a]=vp[a]/s2;
    for(int a=0;a<J;++a) for(int b=0;b<J;++b) Vc[(size_t)a*J+b]-=vp[a]*vp[b]/s2;
  }
  // NOTE: count measurement is skipped on structural-zero particles. Therefore the
  // fixed Kc sequence above is exact only if every count latent is conditioned on.
  // We consequently maintain particle-specific diagonal/full V below; Kc/sc are
  // retained only as harmless precomputed starting diagnostics and not used.

  const double tinyu=std::max(tiny,std::numeric_limits<double>::min());
  const int nt=nth(n_threads,n); std::vector<double> out(n,R_NegInf); std::vector<std::thread> th; th.reserve(nt);

  auto worker=[&](int lo,int hi){
    std::vector<double> mc(J), V((size_t)J*J), vpsi(J), probs(d), mass(d);
    for(int i=lo;i<hi;++i){ double lse=R_NegInf;
      for(int q=0;q<Q;++q){
        std::fill(mc.begin(),mc.end(),0.0); std::fill(V.begin(),V.end(),0.0); for(int a=0;a<J;++a)V[(size_t)a*J+a]=1.0;
        double mb=0.0, vb=1.0, lw=0.0; bool ok=true;
        for(int j=0;j<k;++j){
          const bool iszi=zip[j], is0=iszi && yzero[(size_t)i*k+j];
          double pactive=1.0, sz=1.0, kz=0.0, muz=0.0;
          if(iszi){
            if(L==0) pactive=std::min(1.0-1e-14,std::max(1e-14,pb[j][i]));
            else {
              const double rr=rho[zipos[j]], Dz=std::max(1e-12,1.0-rr*rr);
              muz=rr*mb; double sz2=Dz+rr*rr*vb; sz=std::sqrt(std::max(1e-12,sz2)); kz=vb*rr/sz2;
              pactive=std::min(1.0-1e-14,std::max(1e-14,Phi((qpb[j][i]-muz)/sz)));
            }
          }

          // Current count predictive distribution from particle-specific state.
          double muc=0.0; for(int a=0;a<J;++a)muc+=psi[j+(size_t)k*a]*mc[a];
          for(int a=0;a<J;++a){ double z=0; for(int b=0;b<J;++b)z+=V[(size_t)a*J+b]*psi[j+(size_t)k*b]; vpsi[a]=z; }
          double sc2=Dc[j]; for(int a=0;a<J;++a)sc2+=psi[j+(size_t)k*a]*vpsi[a]; sc2=std::max(1e-12,sc2); double scc=std::sqrt(sc2);

          auto count_norm=[&](const std::vector<double>& Imat){ double z=0; for(int m=0;m<d;++m){ double Fl=(m==0)?0.0:Phi((edge[m]-muc)/scc); double Fu=(m==d-1)?1.0:Phi((edge[m+1]-muc)/scc); probs[m]=std::max(0.0,Fu-Fl); mass[m]=d*std::max(Imat[i+(size_t)n*m],0.0)*probs[m]; z+=mass[m]; } return z; };

          bool active = true;
          double normc = 0.0;
          if(is0){
            normc=count_norm(I0[j]); const double wa=1.0-pactive, wb=pactive*normc, total=wa+wb;
            if(!(total>tinyu) || !std::isfinite(total)){ok=false;break;} lw+=std::log(total);
            double uz=clampu(U_zero(q,j)), branch_target=uz*total; active=(branch_target>wa);
            if(!active){
              // structural zero: condition zero latent to inactive region (> threshold), no count update
              if(L==1){ double within0=clampu(branch_target/std::max(wa,tinyu)); double p=pactive+within0*(1.0-pactive);
                double z=qnorm_safe(p,muz,sz), inn=z-muz; const double rr=rho[zipos[j]]; mb+=kz*inn; vb=std::max(1e-12,vb-kz*(rr*vb)); }
              continue;
            }
            // Re-map the same mixture uniform to U(0,1) within the active branch.
            // This is the Genz conditional transform and avoids branch-selection bias.
            if(L==1){ double within1=clampu((branch_target-wa)/std::max(wb,tinyu)); double p=within1*pactive;
              double z=qnorm_safe(p,muz,sz), inn=z-muz; const double rr=rho[zipos[j]]; mb+=kz*inn; vb=std::max(1e-12,vb-kz*(rr*vb)); }
           // active + sampling zero; normc already computed
          } else {
            normc=count_norm(Iy[j]); double total=pactive*normc;
            if(!(total>tinyu)||!std::isfinite(total)){ok=false;break;} lw+=std::log(total);
           }

          // Positive ZIPM observation is necessarily in the active zero branch.
          // For y=0 the active-branch update was already performed above using
          // the branch-remapped uniform.
          if(iszi && L==1 && !is0){ double p=clampu(U_zero(q,j))*pactive; double z=qnorm_safe(p,muz,sz), inn=z-muz; const double rr=rho[zipos[j]]; mb+=kz*inn; vb=std::max(1e-12,vb-kz*(rr*vb)); }

          // Sample tilted count latent interval/position and update count-factor state.
          if(!(normc>tinyu)){ok=false;break;} if(!is0){} // normc already in lw for positive; for zero active total included normc.
          double target=clampu(U_count(q,j))*normc, acc=0.0; int mm=d-1;
          for(int m=0;m<d;++m){ if(target<=acc+mass[m]){mm=m;break;} acc+=mass[m]; }
          double within=clampu((target-acc)/std::max(mass[mm],tinyu));
          double Fl=(mm==0)?0.0:Phi((edge[mm]-muc)/scc), Fu=(mm==d-1)?1.0:Phi((edge[mm+1]-muc)/scc);
          double x=qnorm_safe(Fl+within*(Fu-Fl),muc,scc), inn=x-muc;
          for(int a=0;a<J;++a) mc[a]+=vpsi[a]/sc2*inn;
          for(int a=0;a<J;++a)for(int b=0;b<J;++b)V[(size_t)a*J+b]-=vpsi[a]*vpsi[b]/sc2;
        }
        if(ok && std::isfinite(lw)) lse=logadd(lse,lw);
      }
      out[i]=lse-std::log((double)Q);
    }
  };
  int base=n/nt,rem=n%nt,st=0; for(int t=0;t<nt;++t){int len=base+(t<rem),lo=st,hi=st+len;st=hi;th.emplace_back(worker,lo,hi);} for(auto &x:th)x.join();
  NumericVector ans(n); for(int i=0;i<n;++i)ans[i]=out[i]; return ans;
}

// Backward-compatible PM-only wrapper using the same clean kernel logic is kept
// in R rather than duplicated here.  This file intentionally exports one final
// likelihood kernel only.
