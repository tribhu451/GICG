//////////////////////////////////////////////////////////////////
//   Impact parameter Dependent Nucleon Nucleon Interaction     //
//////////////////////////////////////////////////////////////////

#pragma once
#include <fstream>
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <TRandom3.h>
#include <TF1.h>
#include "TMath.h"
#include "input_data.h"
#include "mc_glau.h"
#include "cell.h"
#include "grid.h"
#include <gsl/gsl_sf_erf.h>


#define PI  3.1415927

using std::cout;
using std::endl;

class idnni{
  
 public :
  idnni(InputData *inparams_,grid* ,  mc_glau* );
  ~idnni();
  void perform();
  
  double get_total_deposition_over_all_cells_for_current_event();
  void   update_contribution_on_cells_over_all_events() ; 
  void   write_event_averaged_profile_to_file(int nEvents, int flag_to_generate_music_boost_invariant_file);
  void   get_eccentricities_and_participant_plane_angles(int Norder, double* eps, double* psi);

  // gaussian smearing
  void   do_gaussian_smearing_for_the_current_event(double sigma_perp);
  void   rotation_by_second_order_event_plane_followed_by_gaussian_smearing_for_the_current_event(double sigma_perp, double psi2);
  void   update_contribution_on_cells_over_all_events_with_gaussian_smearing();
  void   get_eccentricities_and_participant_plane_angles_after_gaussian_smearing(int Norder, double* eps, double* psi); 
  void   write_event_averaged_profile_to_file_after_gaussian_smearing(int nEvents, int flag_to_generate_music_boost_invariant_file);
  double get_total_deposition_over_all_cells_after_gaussian_smearing_for_current_event();

 private :
  
  mc_glau* mc;
  grid* arena ; 
  InputData *inparams;
  // 2D grid, packed in 1D array
  
  // grid info
  float x_max;
  float y_max;
  int nx;
  int ny;
  float dx;
  float dy;
  void get_nearest_cell_index(double xx, double yy, int& ix, int& iy ) ; 
  void get_cm_of_the_distribution(double& xx, double& yy) ;
  void get_cm_of_the_distribution_after_gaussian_smearing(double& xx, double& yy); 

  
  // XA, YA, ZA -> (x,y,z) coordinate of nucleus A.
  // XB, YB, ZB -> (x,y,z) coordinate of nucleus B.
  double XA[300];
  double YA[300];
  double ZA[300];
  double XB[300];
  double YB[300];
  double ZB[300];
  int pflagA[300]; // proton flag
  int pflagB[300] ; 
  
  double max(double a, double b){
    if(a>b){
      return a;
    }
    else{
      return b;
    }
  }
  
};
