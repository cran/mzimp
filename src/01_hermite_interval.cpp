#include <Rcpp.h>
#include <vector>
#include <cmath>
#include <algorithm>
#include <limits>
#include <thread>
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
inline void parallel_chunks(int nwork, int requested, Func fn) {
  const int nt = effective_threads(requested, nwork);
  if (nt <= 1) { fn(0, nwork, 0); return; }
  const int chunk = (nwork + nt - 1) / nt;
  std::vector<std::thread> workers;
  workers.reserve(nt);
  for (int t = 0; t < nt; ++t) {
    const int begin = t * chunk;
    const int end = std::min(nwork, begin + chunk);
    if (begin >= end) break;
    workers.emplace_back([&, begin, end, t]() { fn(begin, end, t); });
  }
  for (auto& th : workers) th.join();
}
}

namespace {
inline double normal_pdf_std(double x) {
  static const double INV_SQRT_2PI = 0.39894228040143267793994605993438;
  return INV_SQRT_2PI * std::exp(-0.5*x*x);
}
inline double normal_cdf_std(double x) {
  static const double INV_SQRT2 = 0.70710678118654752440084436210485;
  return 0.5 * std::erfc(-x * INV_SQRT2);
}
inline double normal_sf_std(double x) {
  static const double INV_SQRT2 = 0.70710678118654752440084436210485;
  return 0.5 * std::erfc(x * INV_SQRT2);
}
inline double log_normal_cdf_std(double x) {
  if (x < -10.0) {
    const double z=-x, z2=z*z;
    double corr=1.0 - 1.0/z2 + 3.0/(z2*z2) - 15.0/(z2*z2*z2);
    corr=std::max(corr,1e-12);
    return -0.5*x*x - std::log(z) - 0.91893853320467274178032973640562 + std::log(corr);
  }
  const double p=normal_cdf_std(x);
  return p>0.0 ? std::log(p) : -std::numeric_limits<double>::infinity();
}
inline double log_normal_sf_std(double x) { return log_normal_cdf_std(-x); }
inline double clamp01(double x) {
  if (x < 0.0) return 0.0;
  if (x > 1.0) return 1.0;
  return x;
}
inline double logspace_sub_cpp(double logx, double logy) {
  if (!std::isfinite(logx)) return R_NegInf;
  if (!std::isfinite(logy)) return logx;
  if (!(logx > logy)) return R_NegInf;
  const double d = logy - logx;
  return logx + std::log(-std::expm1(d));
}
inline double log_normal_interval(double lo, double up) {
  if (R_IsNA(lo) || R_IsNA(up) || lo >= up) return R_NegInf;
  if (!std::isfinite(lo) && lo < 0.0)
    return log_normal_cdf_std(up);
  if (!std::isfinite(up) && up > 0.0)
    return log_normal_sf_std(lo);
  if (up <= 0.0) {
    return logspace_sub_cpp(
      log_normal_cdf_std(up),
      log_normal_cdf_std(lo)
    );
  }
  if (lo >= 0.0) {
    return logspace_sub_cpp(
      log_normal_sf_std(lo),
      log_normal_sf_std(up)
    );
  }
  const double p = normal_cdf_std(up) - normal_cdf_std(lo);
  return std::log(std::max(p, std::numeric_limits<double>::min()));
}
inline double signed_log_sum_to_value(double x, double dx, double y, double dy) {
  // phi(x)*dx - phi(y)*dy, evaluated in signed-log form.
  double l1 = R_NegInf, l2 = R_NegInf;
  int s1 = 0, s2 = 0;
  if (std::isfinite(x) && dx != 0.0) {
    l1 = std::log(normal_pdf_std(x)) + std::log(std::abs(dx));
    s1 = (dx > 0.0) ? 1 : -1;
  }
  if (std::isfinite(y) && dy != 0.0) {
    l2 = std::log(normal_pdf_std(y)) + std::log(std::abs(dy));
    s2 = (dy > 0.0) ? -1 : 1; // subtract lower-bound term
  }
  if (s1 == 0 && s2 == 0) return 0.0;
  if (s2 == 0) return s1 * std::exp(l1);
  if (s1 == 0) return s2 * std::exp(l2);
  if (s1 == s2) {
    const double a = std::max(l1,l2), b = std::min(l1,l2);
    return s1 * std::exp(a) * (1.0 + std::exp(b-a));
  }
  if (l1 == l2) return 0.0;
  if (l1 > l2) return s1 * std::exp(l1) * (-std::expm1(l2-l1));
  return s2 * std::exp(l2) * (-std::expm1(l1-l2));
}
inline double logsum2_cpp(double a, double b) {
  const double m = std::max(a,b);
  if (!std::isfinite(m)) return m;
  return m + std::log(std::exp(a-m)+std::exp(b-m));
}
}

// -----------------------------------------------------------------------------
// Kernel 1: all checkerboard conditional cell kernels for j=1,...,k.
// This replaces compute_H_kernel()'s R loop over margins.
// -----------------------------------------------------------------------------
// [[Rcpp::export]]
Rcpp::List hermite_interval_kernel_cpp(
    const Rcpp::NumericMatrix& Psi,
    const Rcpp::NumericVector& herm_x,
    const Rcpp::NumericVector& t_inner,
    const int d0,
    int n_threads = 1) {

  const int k = Psi.nrow();
  const int J = Psi.ncol();
  const int Q = herm_x.size();
  if (d0 < 2 || t_inner.size() != d0-1) stop("hermite_interval_kernel_cpp: invalid d0/t_inner.");
  if (Q < 1 || J < 1 || k < 1) stop("hermite_interval_kernel_cpp: invalid dimensions.");

  NumericVector rvec(k), svec(k);
  List OmegaT_list(k), logOmegaT_list(k), dlogOmegaT_dr_list(k), dOmegaT_dr_list(k), dOmegaT_dx_list(k);
  std::vector<NumericMatrix> OmM, lOmM, dlOmM, dOmM, dOmXM;
  OmM.reserve(k); lOmM.reserve(k); dlOmM.reserve(k); dOmM.reserve(k); dOmXM.reserve(k);
  std::vector<double*> OmP(k), lOmP(k), dlOmP(k), dOmP(k), dOmXP(k);

  for (int j=0;j<k;++j) {
    double rr2=0.0;
    for(int a=0;a<J;++a) rr2 += Psi(j,a)*Psi(j,a);
    const double rr=std::sqrt(rr2), ss=std::sqrt(std::max(1e-14,1.0-rr2));
    if(!std::isfinite(rr)||!std::isfinite(ss)||ss<=0.0||rr>=1.0) stop("hermite_interval_kernel_cpp: invalid loading row.");
    rvec[j]=rr; svec[j]=ss;
    NumericMatrix Om(d0,Q), lOm(d0,Q), dlOm(d0,Q), dOm(d0,Q), dOmX(d0,Q);
    OmegaT_list[j]=Om; logOmegaT_list[j]=lOm; dlogOmegaT_dr_list[j]=dlOm; dOmegaT_dr_list[j]=dOm; dOmegaT_dx_list[j]=dOmX;
    OmM.push_back(Om); lOmM.push_back(lOm); dlOmM.push_back(dlOm); dOmM.push_back(dOm); dOmXM.push_back(dOmX);
    OmP[j]=REAL(Om); lOmP[j]=REAL(lOm); dlOmP[j]=REAL(dlOm); dOmP[j]=REAL(dOm); dOmXP[j]=REAL(dOmX);
  }

  const int nt=zipm_parallel::effective_threads(n_threads,k);
  zipm_parallel::parallel_chunks(k, nt, [&](int j_begin, int j_end, int /*tid*/) {
  for(int j=j_begin;j<j_end;++j) {
    const double r=rvec[j], s=svec[j], s3=s*s*s;
    double* Om=OmP[j]; double* lOm=lOmP[j]; double* dlOm=dlOmP[j]; double* dOm=dOmP[j]; double* dOmX=dOmXP[j];
    for(int q=0;q<Q;++q) {
      const double x=herm_x[q];
      for(int m=0;m<d0;++m) {
        double lo=R_NegInf, up=R_PosInf, dlo=0.0, dup=0.0;
        double dlo_x=0.0, dup_x=0.0;
        if(m>0) {
          const double num=t_inner[m-1]-r*x;
          lo=num/s;
          dlo=-x/s + num*r/s3;
          dlo_x=-r/s;
        }
        if(m<d0-1) {
          const double num=t_inner[m]-r*x;
          up=num/s;
          dup=-x/s + num*r/s3;
          dup_x=-r/s;
        }
        const double lom=log_normal_interval(lo,up);
        const double om=std::isfinite(lom) ? std::exp(lom) : 0.0;
        const double dom=signed_log_sum_to_value(up,dup,lo,dlo);
        const double domx=(std::isfinite(up)?normal_pdf_std(up)*dup_x:0.0)
                         -(std::isfinite(lo)?normal_pdf_std(lo)*dlo_x:0.0);
        const int id=m+d0*q;
        Om[id]=om; lOm[id]=lom; dOm[id]=std::isfinite(dom)?dom:0.0; dOmX[id]=std::isfinite(domx)?domx:0.0;
        dlOm[id]=(om>0.0 && std::isfinite(dom)) ? dom/om : 0.0;
      }
    }
  }
  });

  return List::create(
    _["Psi"]=Psi, _["r"]=rvec, _["s"]=svec,
    _["OmegaT_list"]=OmegaT_list,
    _["logOmegaT_list"]=logOmegaT_list,
    _["dlogOmegaT_dr_list"]=dlogOmegaT_dr_list,
    _["dOmegaT_dr_list"]=dOmegaT_dr_list,
    _["dOmegaT_dx_list"]=dOmegaT_dx_list
  );
}

// -----------------------------------------------------------------------------
// Kernel 2: count cache for joint L>=1 Hermite path.
// Replaces make_hermite_count_cache_block() R loop over k.
// -----------------------------------------------------------------------------
// [[Rcpp::export]]
Rcpp::List hermite_count_cache_cpp(
    const Rcpp::List& Iy_list,
    const Rcpp::List& log_sy_list,
    const Rcpp::List& I0_list,
    const Rcpp::List& log_s0_list,
    const Rcpp::List& OmegaT_list,
    const Rcpp::List& dOmegaT_dr_list,
    const Rcpp::List& dOmegaT_dx_list,
    const Rcpp::LogicalVector& zero_inflation,
    const Rcpp::IntegerVector& idx1,
    const double d0,
    const bool need_grad,
    const double tiny=1e-300,
    int n_threads=1) {

  const int k=Iy_list.size(), nb=idx1.size();
  if(zero_inflation.size()!=k) stop("hermite_count_cache_cpp: zero_inflation mismatch.");
  std::vector<int> ridx(nb); for(int i=0;i<nb;++i) ridx[i]=idx1[i]-1;
  List logY(k), dlogY(k), dlogYdx(k), log0(k), dlog0(k), dlog0dx(k);

  struct MRef { int nr,dd,Q; const double *Iy,*sy,*I0,*s0,*Om,*dOm,*dOmX; double *LY,*DY,*DYX,*L0,*D0,*D0X; bool zi; };
  std::vector<MRef> M(k);
  std::vector<NumericMatrix> keepMat; keepMat.reserve(k*10);
  std::vector<NumericVector> keepVec; keepVec.reserve(k*2);

  for(int j=0;j<k;++j) {
    NumericMatrix Iy=Iy_list[j], Om=OmegaT_list[j]; NumericVector sy=log_sy_list[j];
    if(Iy.ncol()!=Om.nrow()||sy.size()!=Iy.nrow()) stop("hermite_count_cache_cpp: dimension mismatch.");
    const int Q=Om.ncol(); NumericMatrix LY(nb,Q); logY[j]=LY;
    NumericMatrix DY, DYX;
    const double *dOmP=nullptr,*dOmXP=nullptr; double *DYP=nullptr,*DYXP=nullptr;
    if(need_grad) {
      NumericMatrix dOm=dOmegaT_dr_list[j], dOmX=dOmegaT_dx_list[j];
      DY=NumericMatrix(nb,Q); DYX=NumericMatrix(nb,Q); dlogY[j]=DY; dlogYdx[j]=DYX;
      dOmP=REAL(dOm); dOmXP=REAL(dOmX); DYP=REAL(DY); DYXP=REAL(DYX); keepMat.push_back(dOm); keepMat.push_back(dOmX); keepMat.push_back(DY); keepMat.push_back(DYX);
    } else { dlogY[j]=R_NilValue; dlogYdx[j]=R_NilValue; }
    const bool zi=zero_inflation[j]; const double *I0P=nullptr,*s0P=nullptr; double *L0P=nullptr,*D0P=nullptr,*D0XP=nullptr;
    if(zi) {
      NumericMatrix I0=I0_list[j]; NumericVector s0=log_s0_list[j]; NumericMatrix L0(nb,Q); log0[j]=L0;
      I0P=REAL(I0); s0P=REAL(s0); L0P=REAL(L0); keepMat.push_back(I0); keepVec.push_back(s0); keepMat.push_back(L0);
      if(need_grad) { NumericMatrix D0(nb,Q),D0X(nb,Q); dlog0[j]=D0; dlog0dx[j]=D0X; D0P=REAL(D0); D0XP=REAL(D0X); keepMat.push_back(D0); keepMat.push_back(D0X); }
      else { dlog0[j]=R_NilValue; dlog0dx[j]=R_NilValue; }
    } else { log0[j]=R_NilValue; dlog0[j]=R_NilValue; dlog0dx[j]=R_NilValue; }
    keepMat.push_back(Iy); keepVec.push_back(sy); keepMat.push_back(Om); keepMat.push_back(LY);
    M[j]={Iy.nrow(),Iy.ncol(),Q,REAL(Iy),REAL(sy),I0P,s0P,REAL(Om),dOmP,dOmXP,REAL(LY),DYP,DYXP,L0P,D0P,D0XP,zi};
  }
  const double logtiny=std::log(tiny); const int nt=zipm_parallel::effective_threads(n_threads,k);
  zipm_parallel::parallel_chunks(k, nt, [&](int j_begin, int j_end, int /*tid*/) {
  for(int j=j_begin;j<j_end;++j) {
    const MRef& a=M[j];
    for(int i=0;i<nb;++i) {
      const int rr=ridx[i];
      for(int q=0;q<a.Q;++q) {
        double py=0,dy=0,dyx=0,p0=0,d0r=0,d0x=0;
        for(int m=0;m<a.dd;++m) {
          const double om=a.Om[m+a.dd*q];
          py += a.Iy[rr+a.nr*m]*om;
          if(a.zi) p0 += a.I0[rr+a.nr*m]*om;
          if(need_grad) {
            const double dor=a.dOm[m+a.dd*q], dox=a.dOmX[m+a.dd*q];
            dy += a.Iy[rr+a.nr*m]*dor; dyx += a.Iy[rr+a.nr*m]*dox;
            if(a.zi) { d0r += a.I0[rr+a.nr*m]*dor; d0x += a.I0[rr+a.nr*m]*dox; }
          }
        }
        py*=d0; dy*=d0; dyx*=d0; const int id=i+nb*q;
        if(std::isfinite(py)&&py>0) { a.LY[id]=std::log(py)+a.sy[rr]; if(need_grad){a.DY[id]=dy/py; a.DYX[id]=dyx/py;} }
        else { a.LY[id]=logtiny+a.sy[rr]; if(need_grad){a.DY[id]=0; a.DYX[id]=0;} }
        if(a.zi) {
          p0*=d0; d0r*=d0; d0x*=d0;
          if(std::isfinite(p0)&&p0>0) { a.L0[id]=std::log(p0)+a.s0[rr]; if(need_grad){a.D0[id]=d0r/p0; a.D0X[id]=d0x/p0;} }
          else { a.L0[id]=logtiny+a.s0[rr]; if(need_grad){a.D0[id]=0; a.D0X[id]=0;} }
        }
      }
    }
  }
  });
  return List::create(_["logY"]=logY,_["dlogY"]=dlogY,_["dlogY_dx"]=dlogYdx,_["log0"]=log0,_["dlog0"]=dlog0,_["dlog0_dx"]=dlog0dx);
}

// -----------------------------------------------------------------------------
// Kernel 3: complete L=0 block preparation + Hermite projection.
// It fuses I*Omega, independent-ZI mixture, radial derivative and c0/c1/c2.
// -----------------------------------------------------------------------------
// [[Rcpp::export]]
Rcpp::List hermite_grid_L0_cpp(
    const Rcpp::List& Iy_list,
    const Rcpp::List& log_sy_list,
    const Rcpp::List& I0_list,
    const Rcpp::List& log_s0_list,
    const Rcpp::List& OmegaT_list,
    const Rcpp::List& dOmegaT_dr_list,
    const Rcpp::List& pb_list,
    const Rcpp::List& y0_list,
    const Rcpp::LogicalVector& zero_inflation,
    const Rcpp::IntegerVector& idx1,
    const Rcpp::NumericMatrix& Wproj,
    const double d0,
    const bool need_grad,
    const double tiny=1e-300,
    int n_threads=1) {

  const int k=Iy_list.size(), nb=idx1.size(), Q=Wproj.nrow();
  if(Wproj.ncol()<3||zero_inflation.size()!=k) stop("hermite_grid_L0_cpp: invalid dimensions.");
  std::vector<int> ridx(nb); for(int i=0;i<nb;++i) ridx[i]=idx1[i]-1;
  NumericMatrix C0(nb,k),C1(nb,k),C2(nb,k),D0(nb,k),D1(nb,k),D2(nb,k);
  List LH(k),dLH(k),FPB(k);

  struct MRef {int nr,dd; const double *Iy,*sy,*I0,*s0,*Om,*dOm,*pb; const int* y0; double *LH,*dLH; bool zi;};
  std::vector<MRef> M(k); std::vector<NumericMatrix> km; std::vector<NumericVector> kv; std::vector<LogicalVector> kl;
  km.reserve(k*7); kv.reserve(k*3); kl.reserve(k);
  for(int j=0;j<k;++j){
    NumericMatrix Iy=Iy_list[j],Om=OmegaT_list[j]; NumericVector sy=log_sy_list[j]; NumericMatrix L(nb,Q); LH[j]=L;
    if(Iy.ncol()!=Om.nrow()||Om.ncol()!=Q) stop("hermite_grid_L0_cpp: Iy/Omega mismatch.");
    const double *dOmP=nullptr; double *dLP=nullptr; if(need_grad){NumericMatrix dOm=dOmegaT_dr_list[j],DL(nb,Q); dLH[j]=DL; dOmP=REAL(dOm); dLP=REAL(DL); km.push_back(dOm); km.push_back(DL);} else dLH[j]=R_NilValue;
    const bool zi=zero_inflation[j]; const double *I0P=nullptr,*s0P=nullptr,*pbP=nullptr; const int* y0P=nullptr;
    if(zi){ NumericMatrix I0=I0_list[j]; NumericVector s0=log_s0_list[j],pb=pb_list[j]; LogicalVector y0=y0_list[j]; I0P=REAL(I0);s0P=REAL(s0);pbP=REAL(pb);y0P=LOGICAL(y0);km.push_back(I0);kv.push_back(s0);kv.push_back(pb);kl.push_back(y0); }
    FPB[j]=R_NilValue; km.push_back(Iy);km.push_back(Om);km.push_back(L);kv.push_back(sy);
    M[j]={Iy.nrow(),Iy.ncol(),REAL(Iy),REAL(sy),I0P,s0P,REAL(Om),dOmP,pbP,y0P,REAL(L),dLP,zi};
  }
  const double *WP=REAL(Wproj); double *c0p=REAL(C0),*c1p=REAL(C1),*c2p=REAL(C2),*d0p=REAL(D0),*d1p=REAL(D1),*d2p=REAL(D2);
  const double logtiny=std::log(tiny); const int nt=zipm_parallel::effective_threads(n_threads,k);
  zipm_parallel::parallel_chunks(k, nt, [&](int j_begin, int j_end, int /*tid*/) {
  for(int j=j_begin;j<j_end;++j){
    const MRef&a=M[j];
    for(int i=0;i<nb;++i){
      const int rr=ridx[i]; double cc0=0,cc1=0,cc2=0,dd0=0,dd1=0,dd2=0;
      for(int q=0;q<Q;++q){
        double py=0,dy=0; for(int m=0;m<a.dd;++m){py+=a.Iy[rr+a.nr*m]*a.Om[m+a.dd*q]; if(need_grad)dy+=a.Iy[rr+a.nr*m]*a.dOm[m+a.dd*q];}
        py*=d0;dy*=d0; double ly=(std::isfinite(py)&&py>0)?std::log(py)+a.sy[rr]:logtiny+a.sy[rr]; double dly=(need_grad&&py>0&&std::isfinite(py))?dy/py:0.0;
        double f=ly,df=dly;
        if(a.zi){ double pb=std::max(1e-12,std::min(1.0-1e-12,a.pb[rr])); if(!a.y0[rr]) f=ly+std::log(pb); else {double p0=0,dp0=0;for(int m=0;m<a.dd;++m){p0+=a.I0[rr+a.nr*m]*a.Om[m+a.dd*q];if(need_grad)dp0+=a.I0[rr+a.nr*m]*a.dOm[m+a.dd*q];}p0*=d0;dp0*=d0;double l0=(p0>0&&std::isfinite(p0))?std::log(p0)+a.s0[rr]:logtiny+a.s0[rr];double dl0=(need_grad&&p0>0&&std::isfinite(p0))?dp0/p0:0.0;double la=std::log1p(-pb),lb=std::log(pb)+l0;f=logsum2_cpp(la,lb);if(need_grad)df=std::exp(lb-f)*dl0;}}
        a.LH[i+nb*q]=f;if(need_grad)a.dLH[i+nb*q]=df;
        const double w0=WP[q+Q*0],w1=WP[q+Q*1],w2=WP[q+Q*2];cc0+=f*w0;cc1+=f*w1;cc2+=f*w2;if(need_grad){dd0+=df*w0;dd1+=df*w1;dd2+=df*w2;}
      }
      c0p[i+nb*j]=cc0;c1p[i+nb*j]=cc1;c2p[i+nb*j]=cc2;if(need_grad){d0p[i+nb*j]=dd0;d1p[i+nb*j]=dd1;d2p[i+nb*j]=dd2;}
    }
  }
  });
  return List::create(_["C0"]=C0,_["C1"]=C1,_["C2"]=C2,_["D0"]=D0,_["D1"]=D1,_["D2"]=D2,_["LH"]=LH,_["dLH"]=dLH,_["FPB"]=FPB,_["dC1_dPB"]=NumericMatrix(nb,k),_["dC2_dPB"]=NumericMatrix(nb,k));
}


// [[Rcpp::export]]
Rcpp::List hermite_count_parallel_info_cpp() {
  const unsigned int hw = std::thread::hardware_concurrency();
  return Rcpp::List::create(
    Rcpp::Named("backend") = "std::thread",
    Rcpp::Named("std_thread") = true,
    Rcpp::Named("hardware_threads") = (hw > 0 ? (int)hw : NA_INTEGER)
  );
}
