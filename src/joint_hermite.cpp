#include <Rcpp.h>
#include <vector>
#include <cmath>
#include <algorithm>
#include <limits>
#include <thread>
#include <atomic>
// [[Rcpp::plugins(cpp11)]]
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

  // Dynamic observation-block scheduling.  We deliberately create many more
  // blocks than workers so expensive observations do not leave other workers
  // idle.  Each block is still a contiguous observation range, hence the
  // lower bound and its gradient remain exactly additive over blocks.
  const int target_blocks = std::max(nt * 4, 1);
  const int block_size = std::max(1, (nwork + target_blocks - 1) / target_blocks);
  std::atomic<int> next_begin(0);

  std::vector<std::thread> workers;
  workers.reserve(nt);
  for (int t = 0; t < nt; ++t) {
    workers.emplace_back([&, t]() {
      while (true) {
        const int begin = next_begin.fetch_add(block_size, std::memory_order_relaxed);
        if (begin >= nwork) break;
        const int end = std::min(nwork, begin + block_size);
        fn(begin, end, t);
      }
    });
  }
  for (auto& th : workers) th.join();
  return nt;
}
}

inline int id2(int r, int c, int D) { return r * D + c; }
inline double clampd(double x, double lo, double hi) { return std::max(lo, std::min(hi, x)); }
inline double normal_pdf_std(double x) {
  static const double INV_SQRT_2PI = 0.39894228040143267793994605993438;
  return INV_SQRT_2PI * std::exp(-0.5*x*x);
}
inline double normal_cdf_std(double x) {
  static const double INV_SQRT2 = 0.70710678118654752440084436210485;
  return 0.5 * std::erfc(-x * INV_SQRT2);
}

static bool chol_decomp(const std::vector<double>& A, std::vector<double>& L, int D) {
  L.assign(D*D, 0.0);
  for (int i=0;i<D;++i) {
    for (int j=0;j<=i;++j) {
      double s=A[id2(i,j,D)];
      for (int k=0;k<j;++k) s -= L[id2(i,k,D)]*L[id2(j,k,D)];
      if (i==j) {
        if (!(s>1e-12) || !R_finite(s)) return false;
        L[id2(i,j,D)] = std::sqrt(s);
      } else {
        L[id2(i,j,D)] = s / L[id2(j,j,D)];
      }
    }
  }
  return true;
}

static void chol_solve(const std::vector<double>& L, const std::vector<double>& b,
                       std::vector<double>& x, int D) {
  std::vector<double> y(D,0.0);
  for(int i=0;i<D;++i){
    double s=b[i];
    for(int k=0;k<i;++k) s -= L[id2(i,k,D)]*y[k];
    y[i]=s/L[id2(i,i,D)];
  }
  x.assign(D,0.0);
  for(int i=D-1;i>=0;--i){
    double s=y[i];
    for(int k=i+1;k<D;++k) s -= L[id2(k,i,D)]*x[k];
    x[i]=s/L[id2(i,i,D)];
  }
}

static void chol_inverse(const std::vector<double>& L, std::vector<double>& V, int D) {
  V.assign(D*D,0.0);
  std::vector<double> e(D), x;
  for(int c=0;c<D;++c){
    std::fill(e.begin(),e.end(),0.0); e[c]=1.0;
    chol_solve(L,e,x,D);
    for(int r=0;r<D;++r) V[id2(r,c,D)] = x[r];
  }
}

inline double logsum2(double a, double b) {
  double m=std::max(a,b);
  if (!R_finite(m)) return m;
  return m + std::log(std::exp(a-m)+std::exp(b-m));
}

// [[Rcpp::export]]
Rcpp::List joint_hermite_cpp(
    Rcpp::List logY_list,
    Rcpp::List dlogY_dr_list,
    Rcpp::List dlogY_dx_list,
    Rcpp::List log0_list,
    Rcpp::List dlog0_dr_list,
    Rcpp::List dlog0_dx_list,
    Rcpp::List y0_list,
    Rcpp::List qpb_list,
    Rcpp::NumericMatrix Psi,
    Rcpp::NumericMatrix Rho,
    Rcpp::IntegerVector zi_position,
    Rcpp::NumericVector herm_x,
    Rcpp::NumericVector herm_w,
    bool need_grad = true,
    double tiny = 1e-300,
    int n_threads = 1
) {
  int k=Psi.nrow(), J=Psi.ncol(), L=Rho.ncol(), D=J+L;
  int Q=herm_x.size();
  if(k==0 || D==0 || Q==0) stop("Invalid dimensions in joint_hermite_allL_cpp().");
  NumericMatrix m0 = logY_list[0];
  int n=m0.nrow();

  // Loading norms and directions.
  std::vector<double> r(k), U(k*J,0.0);
  for(int j=0;j<k;++j){
    double ss=0; for(int a=0;a<J;++a) ss += Psi(j,a)*Psi(j,a);
    r[j]=std::sqrt(ss);
    if(!(r[j]>1e-10) || r[j]>=1.0) stop("Count loading row has invalid norm in joint Hermite backend.");
    for(int a=0;a<J;++a) U[j*J+a]=Psi(j,a)/r[j];
  }
  int nzi=Rho.nrow();
  std::vector<double> t(nzi), VV(nzi*std::max(1,L),0.0);
  for(int s=0;s<nzi;++s){
    double ss=0; for(int a=0;a<L;++a) ss += Rho(s,a)*Rho(s,a);
    t[s]=std::sqrt(ss);
    if(!(t[s]>1e-10) || t[s]>=1.0) stop("Zero loading row has invalid norm in joint Hermite backend.");
    for(int a=0;a<L;++a) VV[s*L+a]=Rho(s,a)/t[s];
  }

  NumericVector log_value(n);
  NumericMatrix gPsi(k,J);
  NumericMatrix gRho(nzi,L);
  const int nt = zipm_parallel::effective_threads(n_threads, n);
  std::vector< std::vector<double> > gPsiThread(nt, std::vector<double>(k*J,0.0));
  std::vector< std::vector<double> > gRhoThread(nt, std::vector<double>(nzi*std::max(1,L),0.0));
  std::vector<double> maxMassThread(nt,0.0), sumMassThread(nt,0.0);
  std::vector<long long> massCountThread(nt,0);

  // Cache matrix/list handles.
  std::vector<NumericMatrix> LY, DYR, DYX;
  std::vector<NumericMatrix> L0, D0R, D0X;
  std::vector<LogicalVector> Y0;
  LY.reserve(k); DYR.reserve(k); DYX.reserve(k); L0.reserve(k); D0R.reserve(k); D0X.reserve(k); Y0.reserve(k);
  for(int j=0;j<k;++j){
    LY.push_back(as<NumericMatrix>(logY_list[j]));
    if(need_grad){ DYR.push_back(as<NumericMatrix>(dlogY_dr_list[j])); DYX.push_back(as<NumericMatrix>(dlogY_dx_list[j])); }
    else { DYR.push_back(NumericMatrix(0,0)); DYX.push_back(NumericMatrix(0,0)); }
    if(zi_position[j]>0){
      L0.push_back(as<NumericMatrix>(log0_list[j]));
      if(need_grad){ D0R.push_back(as<NumericMatrix>(dlog0_dr_list[j])); D0X.push_back(as<NumericMatrix>(dlog0_dx_list[j])); }
      else { D0R.push_back(NumericMatrix(0,0)); D0X.push_back(NumericMatrix(0,0)); }
      Y0.push_back(as<LogicalVector>(y0_list[j]));
    } else {
      L0.push_back(NumericMatrix(0,0)); D0R.push_back(NumericMatrix(0,0)); D0X.push_back(NumericMatrix(0,0)); Y0.push_back(LogicalVector(0));
    }
  }
  std::vector<NumericVector> QPB(k);
  for(int j=0;j<k;++j){ if(zi_position[j]>0) QPB[j]=as<NumericVector>(qpb_list[j]); }

  // Per-observation work. Each thread owns its gradient/diagnostic buffers.
  zipm_parallel::parallel_chunks(n, nt, [&](int i_begin, int i_end, int tid){
    std::vector<double>& gPsiLocal = gPsiThread[tid];
    std::vector<double>& gRhoLocal = gRhoThread[tid];
    double local_max_mass_err=0.0, local_sum_mass_err=0.0;
    long long local_mass_count=0;
  for(int ii=i_begin; ii<i_end; ++ii){
    // coefficient order: c00,c10,c01,c20,c02,c11
    std::vector<double> C(k*6,0.0), Cr(k*6,0.0), Ct(k*6,0.0);
    std::vector<double> avec(D,0.0), B(D*D,0.0);
    double acon=0.0;

    // First pass: standard-normal 1D/2D Hermite projection.
    for(int j=0;j<k;++j){
      int zs=zi_position[j]-1;
      double cj[6]={0,0,0,0,0,0}, crj[6]={0,0,0,0,0,0}, ctj[6]={0,0,0,0,0,0};
      if(zs<0){
        for(int a=0;a<Q;++a){
          double x=herm_x[a], w=herm_w[a], g=LY[j](ii,a);
          double gr = need_grad ? DYR[j](ii,a) : 0.0;
          double bas[6]={1.0,x,0.0,0.5*(x*x-1.0),0.0,0.0};
          for(int h=0;h<6;++h){ cj[h]+=w*g*bas[h]; if(need_grad) crj[h]+=w*gr*bas[h]; }
        }
      } else {
        double tj=t[zs], sig=std::sqrt(std::max(1e-14,1.0-tj*tj)), qpb=QPB[j][ii];
        bool iszero=Y0[j][ii];
        if(!iszero){
          // Positive ZIP cell is exactly separable:
          // g(x,s)=log f_y(x)+log p_B(s).  Under the standard-normal
          // projection measure all mixed Hermite coefficients vanish.
          // Hence the QxQ projection is exactly reducible to two Q loops.
          double ex0=0.0, ex1=0.0, ex2=0.0;
          double er0=0.0, er1=0.0, er2=0.0;
          double es0=0.0, es1=0.0, es2=0.0;
          double et0=0.0, et1=0.0, et2=0.0;
          for(int a=0;a<Q;++a){
            double x=herm_x[a], w=herm_w[a];
            double hx2=0.5*(x*x-1.0);
            double gy=LY[j](ii,a);
            ex0+=w*gy; ex1+=w*gy*x; ex2+=w*gy*hx2;
            if(need_grad){
              double gr=DYR[j](ii,a);
              er0+=w*gr; er1+=w*gr*x; er2+=w*gr*hx2;
            }
          }
          for(int b=0;b<Q;++b){
            double ss=herm_x[b], w=herm_w[b];
            double hs2=0.5*(ss*ss-1.0);
            double aa=(qpb-tj*ss)/sig;
            double pb=clampd(normal_cdf_std(aa),1e-12,1.0-1e-12);
            double gb=std::log(pb);
            es0+=w*gb; es1+=w*gb*ss; es2+=w*gb*hs2;
            if(need_grad){
              double ph=normal_pdf_std(aa);
              double da_dt=-ss/sig + (qpb-tj*ss)*tj/(sig*sig*sig);
              double gt=ph*da_dt/pb;
              et0+=w*gt; et1+=w*gt*ss; et2+=w*gt*hs2;
            }
          }
          cj[0]=ex0+es0; cj[1]=ex1; cj[2]=es1;
          cj[3]=ex2; cj[4]=es2; cj[5]=0.0;
          if(need_grad){
            crj[0]=er0; crj[1]=er1; crj[2]=0.0;
            crj[3]=er2; crj[4]=0.0; crj[5]=0.0;
            ctj[0]=et0; ctj[1]=0.0; ctj[2]=et1;
            ctj[3]=0.0; ctj[4]=et2; ctj[5]=0.0;
          }
        } else {
          // Zero ZIP cell is genuinely coupled in (x,s); keep exact QxQ
          // projection used by v21.
          for(int a=0;a<Q;++a){
            double x=herm_x[a];
            for(int b=0;b<Q;++b){
              double ss=herm_x[b], w=herm_w[a]*herm_w[b];
              double aa=(qpb-tj*ss)/sig;
              double pb=clampd(normal_cdf_std(aa),1e-12,1.0-1e-12);
              double ph=normal_pdf_std(aa);
              double da_dt=-ss/sig + (qpb-tj*ss)*tj/(sig*sig*sig);
              double dpb_dt=ph*da_dt;
              double lp0=L0[j](ii,a);
              double p0=std::exp(std::min(0.0,lp0));
              double lA=std::log1p(-pb), lB=std::log(pb)+lp0, lT=logsum2(lA,lB);
              double g=lT, gr=0.0, gt=0.0;
              if(need_grad){
                double wc=std::exp(lB-lT);
                gr=wc*D0R[j](ii,a);
                double TT=std::max(tiny,std::exp(lT));
                gt=dpb_dt*(p0-1.0)/TT;
              }
              double bas[6]={1.0,x,ss,0.5*(x*x-1.0),0.5*(ss*ss-1.0),x*ss};
              for(int h=0;h<6;++h){
                cj[h]+=w*g*bas[h];
                if(need_grad){ crj[h]+=w*gr*bas[h]; ctj[h]+=w*gt*bas[h]; }
              }
            }
          }
        }
      }
      for(int h=0;h<6;++h){ C[j*6+h]=cj[h]; Cr[j*6+h]=crj[h]; Ct[j*6+h]=ctj[h]; }
      acon += cj[0]-cj[3]-cj[4];
      // linear count block
      for(int a=0;a<J;++a) avec[a] += cj[1]*U[j*J+a];
      // count quadratic block B, exponent is a + d'S + .5 S' B S
      for(int a=0;a<J;++a) for(int b=0;b<J;++b)
        B[id2(a,b,D)] += 2.0*cj[3]*U[j*J+a]*U[j*J+b];
      if(zs>=0){
        for(int a=0;a<L;++a) avec[J+a] += cj[2]*VV[zs*L+a];
        for(int a=0;a<L;++a) for(int b=0;b<L;++b)
          B[id2(J+a,J+b,D)] += 2.0*cj[4]*VV[zs*L+a]*VV[zs*L+b];
        for(int a=0;a<J;++a) for(int b=0;b<L;++b){
          double vv0=cj[5]*U[j*J+a]*VV[zs*L+b];
          B[id2(a,J+b,D)] += vv0; B[id2(J+b,a,D)] += vv0;
        }
      }
    }

    std::vector<double> K(D*D,0.0);
    for(int a=0;a<D;++a) for(int b=0;b<D;++b) K[id2(a,b,D)] = (a==b?1.0:0.0)-B[id2(a,b,D)];
    std::vector<double> chol;
    if(!chol_decomp(K,chol,D)){ log_value[ii]=R_NegInf; continue; }
    std::vector<double> mu, V;
    chol_solve(chol,avec,mu,D); chol_inverse(chol,V,D);
    double logdet=0.0; for(int a=0;a<D;++a) logdet += 2.0*std::log(chol[id2(a,a,D)]);
    double dmu=0.0; for(int a=0;a<D;++a) dmu += avec[a]*mu[a];
    double logbase=acon-0.5*logdet+0.5*dmu;

    double ER=0.0;
    std::vector<double> ERS(D,0.0), ERSS(D*D,0.0);
    std::vector<double> directPsi(k*J,0.0), directRho(nzi*std::max(1,L),0.0);

    // Second pass: first-cumulant residual under tilted Gaussian q.
    for(int j=0;j<k;++j){
      int zs = zi_position[j] - 1;
      bool zip = (zs >= 0);
      std::vector<double> ax(D,0.0), as(D,0.0);
      for(int a=0;a<J;++a) ax[a]=U[j*J+a];
      if(zip) for(int a=0;a<L;++a) as[J+a]=VV[zs*L+a];
      double mx=0.0, ms=0.0;
      for(int a=0;a<D;++a){ mx += ax[a]*mu[a]; if(zip) ms += as[a]*mu[a]; }
      auto qform=[&](const std::vector<double>& aa,const std::vector<double>& bb){
        double z=0.0; for(int a=0;a<D;++a) for(int b=0;b<D;++b) z += aa[a]*V[id2(a,b,D)]*bb[b]; return z; };
      double vx=std::max(1e-10,qform(ax,ax));
      double vs=zip?std::max(1e-10,qform(as,as)):1.0;
      double cxs=zip?qform(ax,as):0.0;
      double det=zip?(vx*vs-cxs*cxs):vx;
      if(!(det>1e-14)){ log_value[ii]=R_NegInf; break; }
      // cov(S,Y)
      std::vector<double> covx(D,0.0), covs(D,0.0);
      for(int a=0;a<D;++a){
        for(int b=0;b<D;++b){ covx[a]+=V[id2(a,b,D)]*ax[b]; if(zip) covs[a]+=V[id2(a,b,D)]*as[b]; }
      }
      // Kloc = cov(S,Y) Sig^{-1}
      std::vector<double> kx(D,0.0), ks(D,0.0);
      double i11,i12=0.0,i22=0.0;
      if(zip){ i11=vs/det; i22=vx/det; i12=-cxs/det; }
      else i11=1.0/vx;
      for(int a=0;a<D;++a){
        kx[a]=covx[a]*i11 + (zip?covs[a]*i12:0.0);
        if(zip) ks[a]=covx[a]*i12+covs[a]*i22;
      }
      // Conditional covariance Vc = V - Kloc Cov(Y,S)
      std::vector<double> Vc=V;
      for(int a=0;a<D;++a) for(int b=0;b<D;++b){
        Vc[id2(a,b,D)] -= kx[a]*covx[b] + (zip?ks[a]*covs[b]:0.0);
      }

      // One-pass normalized fixed-grid expectation.
      // Algebraically identical to the original maxlw -> normw -> moment
      // computation, but each node is visited only once.
      double h0=0.0, hd1=0.0, hd2=0.0, h11=0.0, h12=0.0, h22=0.0;
      double Egr=0.0,Egt=0.0,Egx=0.0,Egs=0.0;
      double Egxd1=0.0,Egxd2=0.0,Egsd1=0.0,Egsd2=0.0;
      double cj0=C[j*6+0], cj1=C[j*6+1], cj2=C[j*6+2],
             cj3=C[j*6+3], cj4=C[j*6+4], cj5=C[j*6+5];
      double tj=zip?t[zs]:0.0;
      double sig=zip?std::sqrt(std::max(1e-14,1.0-tj*tj)):1.0;
      double qpb=zip?QPB[j][ii]:0.0;
      bool iszero=zip?Y0[j][ii]:false;

      double lmax=-std::numeric_limits<double>::infinity(), normw=0.0;
      auto scale_sums=[&](double sc){
        h0*=sc; hd1*=sc; hd2*=sc; h11*=sc; h12*=sc; h22*=sc;
        Egr*=sc; Egt*=sc; Egx*=sc; Egs*=sc;
        Egxd1*=sc; Egxd2*=sc; Egsd1*=sc; Egsd2*=sc;
      };

      for(int a=0;a<Q;++a){
        double x=herm_x[a], dx=x-mx;
        int bmax=zip?Q:1;
        for(int bb=0;bb<bmax;++bb){
          double ss=zip?herm_x[bb]:0.0, ds=zip?(ss-ms):0.0;
          double lr,lwb;
          if(!zip){
            lr=-0.5*std::log(vx)-0.5*dx*dx/vx+0.5*x*x;
            lwb=std::log(std::max(herm_w[a],tiny))+lr;
          } else {
            double quad=(vs*dx*dx-2*cxs*dx*ds+vx*ds*ds)/det;
            lr=-0.5*std::log(det)-0.5*quad+0.5*(x*x+ss*ss);
            lwb=std::log(std::max(herm_w[a]*herm_w[bb],tiny))+lr;
          }

          double ww;
          if(!R_finite(lmax)){
            lmax=lwb; normw=1.0; ww=1.0;
          } else if(lwb<=lmax){
            ww=std::exp(lwb-lmax); normw+=ww;
          } else {
            double sc=std::exp(lmax-lwb);
            scale_sums(sc);
            normw=normw*sc+1.0; lmax=lwb; ww=1.0;
          }

          double g,gx=0.0,gr=0.0,gs=0.0,gt=0.0;
          if(!zip){
            g=LY[j](ii,a);
            if(need_grad){ gx=DYX[j](ii,a); gr=DYR[j](ii,a); }
          } else {
            double aa=(qpb-tj*ss)/sig;
            double pb=clampd(normal_cdf_std(aa),1e-12,1.0-1e-12);
            double ph=normal_pdf_std(aa);
            double dpb_ds=ph*(-tj/sig);
            double da_dt=-ss/sig+(qpb-tj*ss)*tj/(sig*sig*sig);
            double dpb_dt=ph*da_dt;
            if(!iszero){
              g=LY[j](ii,a)+std::log(pb);
              if(need_grad){ gx=DYX[j](ii,a); gr=DYR[j](ii,a); gs=dpb_ds/pb; gt=dpb_dt/pb; }
            } else {
              double lp0=L0[j](ii,a), p0=std::exp(std::min(0.0,lp0));
              double lA=std::log1p(-pb), lB=std::log(pb)+lp0, lT=logsum2(lA,lB);
              g=lT;
              if(need_grad){
                double wc=std::exp(lB-lT), TT=std::max(tiny,std::exp(lT));
                gx=wc*D0X[j](ii,a); gr=wc*D0R[j](ii,a);
                gs=dpb_ds*(p0-1.0)/TT; gt=dpb_dt*(p0-1.0)/TT;
              }
            }
          }
          double P=(cj0-cj3-cj4)+cj1*x+cj2*ss+cj3*x*x+cj4*ss*ss+cj5*x*ss;
          double res=g-P;
          h0+=ww*res; hd1+=ww*res*dx; h11+=ww*res*dx*dx;
          if(zip){ hd2+=ww*res*ds; h12+=ww*res*dx*ds; h22+=ww*res*ds*ds; }
          if(need_grad){
            Egr+=ww*gr; Egx+=ww*gx; Egxd1+=ww*gx*dx;
            if(zip){ Egt+=ww*gt; Egs+=ww*gs; Egxd2+=ww*gx*ds; Egsd1+=ww*gs*dx; Egsd2+=ww*gs*ds; }
          }
        }
      }
      if(!(normw>0.0) || !R_finite(normw)){ log_value[ii]=R_NegInf; break; }

      double invnorm=1.0/normw;
      h0*=invnorm; hd1*=invnorm; hd2*=invnorm; h11*=invnorm; h12*=invnorm; h22*=invnorm;
      Egr*=invnorm; Egt*=invnorm; Egx*=invnorm; Egs*=invnorm;
      Egxd1*=invnorm; Egxd2*=invnorm; Egsd1*=invnorm; Egsd2*=invnorm;

      double logmass=lmax+std::log(normw);
      double mass=(logmass<700.0)?std::exp(logmass):R_PosInf;
      double merr=R_finite(mass)?std::fabs(mass-1.0):R_PosInf;
      local_max_mass_err=std::max(local_max_mass_err,merr);
      local_sum_mass_err+=merr; local_mass_count++;

      ER += h0;
      // ERS contribution.
      std::vector<double> kh(D,0.0);
      for(int a=0;a<D;++a) kh[a]=kx[a]*hd1+(zip?ks[a]*hd2:0.0);
      for(int a=0;a<D;++a) ERS[a] += mu[a]*h0+kh[a];
      // ERSS contribution.
      for(int a=0;a<D;++a) for(int b=0;b<D;++b){
        double kk = kx[a]*h11*kx[b];
        if(zip) kk += kx[a]*h12*ks[b]+ks[a]*h12*kx[b]+ks[a]*h22*ks[b];
        ERSS[id2(a,b,D)] += Vc[id2(a,b,D)]*h0 + mu[a]*mu[b]*h0 + mu[a]*kh[b]+kh[a]*mu[b]+kk;
      }
      if(need_grad){
        // E[g_x S] / E[g_s S] via conditional Gaussian mean.
        std::vector<double> egxS(D,0.0), egsS(D,0.0);
        for(int a=0;a<D;++a){
          egxS[a]=mu[a]*Egx+kx[a]*Egxd1+(zip?ks[a]*Egxd2:0.0);
          if(zip) egsS[a]=mu[a]*Egs+kx[a]*Egsd1+ks[a]*Egsd2;
        }
        // direct loading gradients.
        double dotu=0.0; for(int a=0;a<J;++a) dotu += U[j*J+a]*egxS[a];
        for(int a=0;a<J;++a) directPsi[j*J+a] = Egr*U[j*J+a] + (egxS[a]-U[j*J+a]*dotu)/r[j];
        if(zip){
          double dotv=0.0; for(int a=0;a<L;++a) dotv += VV[zs*L+a]*egsS[J+a];
          for(int a=0;a<L;++a) directRho[zs*L+a] = Egt*VV[zs*L+a] + (egsS[J+a]-VV[zs*L+a]*dotv)/t[zs];
        }
      }
    }
    if(!R_finite(log_value[ii]) && log_value[ii]==R_NegInf) continue;
    log_value[ii]=logbase+ER;
    if(!need_grad) continue;

    // Global residual covariance with S and SS'.
    std::vector<double> gS(D,0.0), H(D*D,0.0);
    for(int a=0;a<D;++a) gS[a]=ERS[a]-ER*mu[a];
    for(int a=0;a<D;++a) for(int b=0;b<D;++b)
      H[id2(a,b,D)] = ERSS[id2(a,b,D)] - ER*(V[id2(a,b,D)]+mu[a]*mu[b]);

    for(int j=0;j<k;++j){
      int zs=zi_position[j]-1; bool zip=zs>=0;
      const double* u=&U[j*J];
      double c10=C[j*6+1], c01=C[j*6+2], c20=C[j*6+3], c02=C[j*6+4], c11=C[j*6+5];
      // contractions.
      double gzu=0.0, uHzu=0.0;
      std::vector<double> Hzu(J,0.0), Hzwv(J,0.0);
      for(int a=0;a<J;++a){
        gzu += gS[a]*u[a];
        for(int b=0;b<J;++b) Hzu[a] += H[id2(a,b,D)]*u[b];
      }
      for(int a=0;a<J;++a) uHzu += u[a]*Hzu[a];
      double gwv=0.0,vHwv=0.0,uHzwv=0.0;
      std::vector<double> Hwv(L,0.0), Hwzu(L,0.0);
      if(zip){
        const double* v=&VV[zs*L];
        for(int a=0;a<L;++a){
          gwv += gS[J+a]*v[a];
          for(int b=0;b<L;++b) Hwv[a] += H[id2(J+a,J+b,D)]*v[b];
          for(int b=0;b<J;++b) Hwzu[a] += H[id2(J+a,b,D)]*u[b];
        }
        for(int a=0;a<L;++a) vHwv += v[a]*Hwv[a];
        for(int a=0;a<J;++a){
          double z=0.0; for(int b=0;b<L;++b) z += H[id2(a,J+b,D)]*v[b]; Hzwv[a]=z; uHzwv += u[a]*z;
        }
      }
      double covrad_r = Cr[j*6+1]*gzu + Cr[j*6+3]*uHzu;
      if(zip) covrad_r += Cr[j*6+2]*gwv + Cr[j*6+4]*vHwv + Cr[j*6+5]*uHzwv;
      std::vector<double> gu(J,0.0);
      for(int a=0;a<J;++a) gu[a]=c10*gS[a]+2.0*c20*Hzu[a]+(zip?c11*Hzwv[a]:0.0);
      double udotgu=0.0; for(int a=0;a<J;++a) udotgu += u[a]*gu[a];
      for(int a=0;a<J;++a){
        double gcov=covrad_r*u[a]+(gu[a]-u[a]*udotgu)/r[j];
        gPsiLocal[j*J+a] += directPsi[j*J+a]+gcov;
      }
      if(zip){
        const double* v=&VV[zs*L];
        double covrad_t = Ct[j*6+1]*gzu + Ct[j*6+2]*gwv + Ct[j*6+3]*uHzu + Ct[j*6+4]*vHwv + Ct[j*6+5]*uHzwv;
        std::vector<double> gv(L,0.0);
        for(int a=0;a<L;++a) gv[a]=c01*gS[J+a]+2.0*c02*Hwv[a]+c11*Hwzu[a];
        double vdot=0.0; for(int a=0;a<L;++a) vdot += v[a]*gv[a];
        for(int a=0;a<L;++a){
          double gcov=covrad_t*v[a]+(gv[a]-v[a]*vdot)/t[zs];
          gRhoLocal[zs*std::max(1,L)+a] += directRho[zs*L+a]+gcov;
        }
      }
    }
  }
    maxMassThread[tid]=std::max(maxMassThread[tid], local_max_mass_err);
    sumMassThread[tid]+=local_sum_mass_err;
    massCountThread[tid]+=local_mass_count;
  });

  double max_mass_err=0.0, sum_mass_err=0.0; long long mass_count=0;
  for(int tid=0; tid<nt; ++tid){
    max_mass_err=std::max(max_mass_err,maxMassThread[tid]);
    sum_mass_err+=sumMassThread[tid]; mass_count+=massCountThread[tid];
    for(int j=0;j<k;++j) for(int a=0;a<J;++a) gPsi(j,a)+=gPsiThread[tid][j*J+a];
    for(int z=0;z<nzi;++z) for(int a=0;a<L;++a) gRho(z,a)+=gRhoThread[tid][z*std::max(1,L)+a];
  }

  return List::create(
    _["log_value"]=log_value,
    _["grad_Psi"]=gPsi,
    _["grad_Rho"]=gRho,
    _["mass_max_error"]=max_mass_err,
    _["mass_mean_error"]=(mass_count>0?sum_mass_err/(double)mass_count:NA_REAL)
  );
}


// [[Rcpp::export]]
Rcpp::List joint_hermite_parallel_info_cpp() {
  const unsigned int hw=std::thread::hardware_concurrency();
  return Rcpp::List::create(
    Rcpp::Named("backend")="std::thread",
    Rcpp::Named("std_thread")=true,
    Rcpp::Named("hardware_threads")=(hw>0?(int)hw:NA_INTEGER),
    Rcpp::Named("max_threads")=(hw>0?(int)hw:1)
  );
}


// -----------------------------------------------------------------------------
// Specialized public wrappers for L = 0 and L = 1.
//
// Design:
//   * L0 has no Bernoulli-factor quadrature/loading and therefore calls the
//     common Hermite core with an empty Rho matrix and all zi_position = 0.
//   * L1 uses the same tested core but validates that Rho has exactly one
//     column.  Keeping separate entry points lets the R layer dispatch once,
//     avoiding L checks in its observation/block loop and leaving room for
//     further L0/L1-specific optimization without changing the API.
// -----------------------------------------------------------------------------

// [[Rcpp::export]]
Rcpp::List joint_hermite_L0_fast_cpp_v22(
    Rcpp::List logY_list,
    Rcpp::List dlogY_dr_list,
    Rcpp::List dlogY_dx_list,
    Rcpp::NumericMatrix Psi,
    Rcpp::NumericVector herm_x,
    Rcpp::NumericVector herm_w,
    bool need_grad = true,
    double tiny = 1e-300,
    int n_threads = 1
) {
  const int k = Psi.nrow();
  if (k <= 0) Rcpp::stop("Psi must have at least one row.");

  Rcpp::List empty_list(k);
  Rcpp::IntegerVector zi_position(k, 0);
  Rcpp::NumericMatrix Rho(0, 0);

  return joint_hermite_cpp(
    logY_list,
    dlogY_dr_list,
    dlogY_dx_list,
    empty_list,
    empty_list,
    empty_list,
    empty_list,
    empty_list,
    Psi,
    Rho,
    zi_position,
    herm_x,
    herm_w,
    need_grad,
    tiny,
    n_threads
  );
}

// [[Rcpp::export]]
Rcpp::List joint_hermite_L1_fast_cpp_v22(
    Rcpp::List logY_list,
    Rcpp::List dlogY_dr_list,
    Rcpp::List dlogY_dx_list,
    Rcpp::List log0_list,
    Rcpp::List dlog0_dr_list,
    Rcpp::List dlog0_dx_list,
    Rcpp::List y0_list,
    Rcpp::List qpb_list,
    Rcpp::NumericMatrix Psi,
    Rcpp::NumericMatrix Rho,
    Rcpp::IntegerVector zi_position,
    Rcpp::NumericVector herm_x,
    Rcpp::NumericVector herm_w,
    bool need_grad = true,
    double tiny = 1e-300,
    int n_threads = 1
) {
  if (Rho.ncol() != 1)
    Rcpp::stop("joint_hermite_L1_fast_cpp_v22 requires Rho with exactly one column (L = 1).");

  return joint_hermite_cpp(
    logY_list,
    dlogY_dr_list,
    dlogY_dx_list,
    log0_list,
    dlog0_dr_list,
    dlog0_dx_list,
    y0_list,
    qpb_list,
    Psi,
    Rho,
    zi_position,
    herm_x,
    herm_w,
    need_grad,
    tiny,
    n_threads
  );
}

// [[Rcpp::export]]
Rcpp::List joint_hermite_L01_fast_info_cpp_v22() {
  const unsigned int hw = std::thread::hardware_concurrency();
  return Rcpp::List::create(
    Rcpp::Named("backend") = "std::thread dynamic observation-block parallel",
    Rcpp::Named("specialized_entry_points") = Rcpp::CharacterVector::create("L0", "L1"),
    Rcpp::Named("shared_core") = "joint_hermite_cpp",
    Rcpp::Named("scheduler") = "dynamic contiguous observation blocks",
    Rcpp::Named("hardware_threads") = (hw > 0 ? (int)hw : NA_INTEGER),
    Rcpp::Named("max_threads") = (hw > 0 ? (int)hw : 1)
  );
}
