// ############################
// # author: Roman Shopa
// # Roman.Shopa[at]ncbj.gov.pl
// ############################
// Functions for the estimation of intensity at any ROR point using 3D kernel
#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

// cumulative distrivution
NumericVector probCDF_cpp(const NumericVector &x,
                          const double & sigma,
                          const double & half_bin_size){
  unsigned int Length = x.length();
  NumericVector Output(Length);
  for(unsigned int i = 0; i < Length; i++) {
    Output[i] = 0.5*(erf((x[i]+half_bin_size)/(sigma*sqrt(2))) -
                     erf((x[i]-half_bin_size)/(sigma*sqrt(2))));
  }
  return Output;
}

// interpolation over discrete high-pass filter profile ('anti-Gauss')
NumericVector antiGaussInterp_cpp(const NumericVector &x,
                                  const NumericVector &int_profile,
                                  const unsigned int &range_index){
  unsigned int x_length = x.length();
  NumericVector Output(x_length);
  for(unsigned int i = 0; i < x_length; i++){
    if(std::abs(x[i]) <= range_index){
      double x_adj = x[i] + range_index; // adjust from symmetric function to indices (0:50)
      if(remainder(x_adj,1)==0) Output[i] = int_profile[int(x_adj)];
      else{
        unsigned int floor_id = int(floor(x_adj));
        unsigned int ceil_id  = int(ceil(x_adj));
        Output[i] = (x_adj-floor_id)*(int_profile[ceil_id]-int_profile[floor_id])/
          (ceil_id-floor_id) + int_profile[floor_id];
      }
    }
  }
  return Output;
}

// choose between three kernels
NumericVector customKernel_cpp(const NumericVector &x,
                               const List          &params, // see below
                               const std::string   &kernel_type){
  if(kernel_type=="cdf") return probCDF_cpp(x, params[0], params[1]);      // 0 - sigma, 1 - half-span
    else if(kernel_type=="gauss") return dnorm(x, 0.0, double(params[0])); // 0 - sigma
    else return antiGaussInterp_cpp(x/params[0], params[1], params[2]);    // 0 - sigma, 1 - dummy profile, 2 - half-span
}

// main function
// [[Rcpp::export]]
NumericVector estimateIntensity_cpp(const NumericVector &x,
                                    const NumericVector &y,
                                    const NumericVector &z,
                                    const NumericVector &ann_point,
                                    const double &cos_Q,  const double &sin_Q,
                                    const double &cos_Fi, const double &sin_Fi,
                                    const List &TOF_params,
                                    const List &z_params,
                                    const List &FBP_params,
                                    const std::string &TOF_filter_type,
                                    const std::string &z_filter_type){
  unsigned int x_length = x.length(); // vectors' length, must be the same for x,y,z
  NumericVector Output(x_length), 
                XFromAnnPt(x_length), YFromAnnPt(x_length), ZFromAnnPt(x_length),
                xr(x_length), Dr(x_length), Dz(x_length), DRL(x_length);
  XFromAnnPt = x - ann_point[0];
  YFromAnnPt = y - ann_point[1];
  ZFromAnnPt = z - ann_point[2];
  xr  = (XFromAnnPt*cos_Fi) + (YFromAnnPt*sin_Fi);
  Dr  = xr*sin_Q + ZFromAnnPt*cos_Q;
  if(sin_Q==0) Dz = 0; else Dz = xr*cos_Q/sin_Q - ZFromAnnPt;
  DRL = -XFromAnnPt*sin_Fi + YFromAnnPt*cos_Fi;
  return customKernel_cpp(Dz, z_params, z_filter_type)* 
         customKernel_cpp(Dr, TOF_params, TOF_filter_type) *
         antiGaussInterp_cpp(DRL/FBP_params[0], FBP_params[1], FBP_params[2]);
}