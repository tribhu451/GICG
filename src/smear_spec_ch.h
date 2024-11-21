#pragma once
#include<iostream>
#include "input_data.h"
#include "mc_glau.h"
#include "cell.h"
#include "grid.h"


class smear_spec_ch {

  public :

    smear_spec_ch(InputData *inparams_,grid* arena_, mc_glau* _mc);
    ~smear_spec_ch();
    void smear_spectator_protons_charge(double , double );
    void update_charge_contribution_on_cells_over_all_events_with_spectator_gaussian_smearing();
    void write_event_averaged_spectator_charge_profile_to_file_after_gaussian_smearing
         (int nEvents, int event_index);
    void write_event_averaged_spectator_charge_profile_to_file_after_gaussian_smearing_with_symmetry
         (int nEvents, int event_index);
    void reset_spectator_charge_contribution_from_all_events_to_zero_on_the_cells();
    void write_eBy_at_000_from_event_averaged_spectator_charge_profile_to_file_assuming_symmetry(int nEvents, int index);

  private :
    mc_glau* mc;
    grid* arena ; 
    InputData *inparams;
    smear_spec_ch* specsm ; 

    int proton_flag_of_the_spectator_in_nucleus_a[500] ; 
    int proton_flag_of_the_spectator_in_nucleus_b[500] ;
    double spectator_x_of_nucleus_a[500] ;
    double spectator_y_of_nucleus_a[500] ;
    double spectator_x_of_nucleus_b[500] ;
    double spectator_y_of_nucleus_b[500] ;

    void get_nearest_cell_index(double xx, double yy, int& ix, int& iy ){
      ix =  ceil( ( xx - (-arena->get_xmax() ) ) / arena->get_dx() ) ; 
      iy =  ceil( ( yy - (-arena->get_ymax() ) ) / arena->get_dy() ) ; 
    }
    double get_beam_rapidity();
    double get_By(double Q, double tau,double etas,double xx,double yy,double xprime,double yprime,double Yb,double emsigma);
    double get_Delta(double tau,double etas,double xx,double yy,double xprime,double yprime,double Yb);
    double get_A(double tau,double etas,double Yb,double Delta,double emsigma);

};
