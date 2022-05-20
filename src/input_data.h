#ifndef SRC_DATA_H
#define SRC_DATA_H

#include<iostream>
#include<string>
using std::string;
 
typedef struct input_data
{

  int mode ;  

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
  int    neta;
  double matter_eta_plateau;
  double matter_eta_fall;
  double matter_tilt_param;
 
  int     baryon_rapidity_profile_type ;
  double  two_component_baryon_deposition_parameter_omega ; 

  double baryon_rapidity_profile_eta_peak;
  double baryon_rapidity_profile_sigma_eta_plus;
  double baryon_rapidity_profile_sigma_eta_minus;

  double iiserbpr_ra_type_2_peak ; 
  double iiserbpr_ra_type_2_etal ; 
  double iiserbpr_ra_type_2_right_fall ; 
 
    
}InputData;

#endif



