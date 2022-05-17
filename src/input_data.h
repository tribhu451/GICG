#ifndef SRC_DATA_H
#define SRC_DATA_H

#include<iostream>
#include<string>
using std::string;
 
typedef struct input_data
{

  double SNN;
  string projectile;
  string target;  

  double bmin,bmax;

  double xhard;
  double npp;

  string root_output_file_name ;

  double xmax ;
  double ymax ;

  int nx ;
  int ny ;

  double gaussian_smearing_sigma ;

  double lower_mult_proxy_cut ;  
  double upper_mult_proxy_cut ; 

  double etamax;
  int neta;
  double matter_eta_plateau;
  double matter_eta_fall;
  double matter_tilt_param;

  double baryon_rapidity_profile_eta_peak;
  double baryon_rapidity_profile_sigma_eta_plus;
  double baryon_rapidity_profile_sigma_eta_minus;
 
    
}InputData;

#endif
