#pragma once
#include<iostream>
#include "input_data.h"
#include "mc_glau.h"
#include "smear_spec_ch.h"
#include "cell.h"
#include "grid.h"


class mc_glau_smear {

  public :

    mc_glau_smear(InputData *inparams_,grid* arena_, mc_glau* _mc);
    ~mc_glau_smear();
    void smear_it(double );
    void update_contribution_on_cells_over_all_events_with_gaussian_smearing();
    void write_event_averaged_profile_to_file_after_gaussian_smearing(int nEvents, int flag_to_generate_music_boost_invariant_file, int event_index);
    void calculate_and_write_EM_field_profile_to_file(int nEvents, int event_index);
    void reset_contribution_from_all_events_to_zero_on_the_cells();
    void write_nucleus_details_of_one_event(int aidx);
  private :
    mc_glau* mc;
    grid* arena ; 
    InputData *inparams;
    smear_spec_ch *specsm ; 

    int    npart, ncoll                   ; 
    double npart_x[500],npart_y[500]      ;
    double npart_x_of_nucleus_a[500],npart_y_of_nucleus_a[500]      ;
    double npart_x_of_nucleus_b[500],npart_y_of_nucleus_b[500]      ;
    double ncoll_x[10000],ncoll_y[10000]  ;
    

    void get_nearest_cell_index(double xx, double yy, int& ix, int& iy ){
      ix =  ceil( ( xx - (-arena->get_xmax() ) ) / arena->get_dx() ) ; 
      iy =  ceil( ( yy - (-arena->get_ymax() ) ) / arena->get_dy() ) ; 
    }



};
