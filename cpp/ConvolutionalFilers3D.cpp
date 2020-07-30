// ############################
// # author: Roman Shopa
// # Roman.Shopa[at]ncbj.gov.pl
// ############################
// Functions for 3D convolutional median and mean average filters
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
#include <cmath>

// calculates median for a vector
double calcMedian(std::vector<double> &in_vector){
  double median;
  size_t size = in_vector.size();
  // sort (stable_sort)
  std::stable_sort(in_vector.begin(),in_vector.end());
  if(size % 2 == 0){
    median = (in_vector[size/2 - 1] + in_vector[size/2]) / 2;
  } else {
    median = in_vector[size/2];
  }
  return median;
}
// calculates mean or median around the voxel with IDs [x_id,y_id,z_id]
double filterVoxel(const unsigned int &x_id,
                   const unsigned int &y_id,
                   const unsigned int &z_id,
                   const arma::cube   &image,
                   const unsigned int &filter_radius,   // filter radius in px
                   const bool         &rounded,         // if true, use baloon not cube
                   const std::string  &filter_type){    // median or mean
  std::vector<double> average_vector = {};
  for(int k = z_id-filter_radius; k < z_id+filter_radius+1; k++){
    for(int i = x_id-filter_radius; i < x_id+filter_radius+1; i++){
      for(int j = y_id-filter_radius; j < y_id+filter_radius+1; j++){
        if(!rounded) average_vector.push_back(image(i,j,k));
        else {
          if(sqrt(pow(int(i-x_id),2) + pow(int(j-y_id),2) + pow(int(k-z_id),2)) <= filter_radius)
            average_vector.push_back(image(i,j,k));
        } 
      }
    }
  }
  if(filter_type=="median") return calcMedian(average_vector); // median filter
  else {
    auto n = average_vector.size(); 
    float average = 0.0f;
    if ( n != 0) {
      average = accumulate(average_vector.begin(), average_vector.end(), 0.0) / n; 
    }
    return average; // mean
  }
}
// filter the whole image  
// [[Rcpp::export]]
arma::cube filterImage_cpp(const arma::cube   &image,
                           const unsigned int &filter_radius,
                           const bool         &rounded,
                           const std::string  &filter_type,
                           const bool         &verbose){    // display progress or not
  // first get image size
  unsigned int xdim = image.n_rows;
  unsigned int ydim = image.n_cols;
  unsigned int zdim = image.n_slices;
  int n_voxels = xdim * ydim * zdim;
  int iter = 0; // iterator for the visualization
  // dummy (empty image)
  arma::cube Output(xdim,ydim,zdim);
  // main loop
  for(unsigned int i = filter_radius; i < (xdim-filter_radius); i++){
    for(unsigned int j = filter_radius; j < (ydim-filter_radius); j++){
      for(unsigned int k = filter_radius; k < (zdim-filter_radius); k++){
        Output(i,j,k) = filterVoxel(i,j,k,image,filter_radius,rounded,filter_type);
        if(verbose){
          iter+=1;
          if(iter % 10000 == 0) std::cout << "Voxels updated: " << iter << "\r";
        }
      }
    }
  }
  return Output;
}
