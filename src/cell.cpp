#include "cell.h"

cell::cell(){
  contri_from_nucleus_a  = 0 ; 
  contri_from_nucleus_b  = 0 ;

  contri_from_nucleus_a_after_gaussian_smearing          = 0 ;
  contri_from_nucleus_b_after_gaussian_smearing          = 0 ;
  contri_from_binary_collisions_after_gaussian_smearing  = 0 ; 

  contri_from_nucleus_a_over_all_events                  = 0 ; 
  contri_from_nucleus_b_over_all_events                  = 0 ;
  contri_from_binary_collisions_over_all_events          = 0 ; 

  charge_contri_from_spectator_nucleus_a_after_gaussian_smearing = 0 ;
  charge_contri_from_spectator_nucleus_b_after_gaussian_smearing = 0;

  charge_contri_from_spectator_nucleus_a_over_all_events = 0 ; 
  charge_contri_from_spectator_nucleus_b_over_all_events = 0 ;


}

void cell::reset_contributions_to_zero(){
  contri_from_nucleus_a  = 0 ; 
  contri_from_nucleus_b  = 0 ;

  contri_from_nucleus_a_after_gaussian_smearing          = 0 ;
  contri_from_nucleus_b_after_gaussian_smearing          = 0 ;
  contri_from_binary_collisions_after_gaussian_smearing  = 0 ; 

  contri_from_nucleus_a_over_all_events                  = 0 ; 
  contri_from_nucleus_b_over_all_events                  = 0 ;
  contri_from_binary_collisions_over_all_events          = 0 ; 
}


void cell::reset_spectator_charge_contributions_to_zero(){
  charge_contri_from_spectator_nucleus_a_after_gaussian_smearing = 0 ;
  charge_contri_from_spectator_nucleus_b_after_gaussian_smearing = 0;

  charge_contri_from_spectator_nucleus_a_over_all_events = 0 ; 
  charge_contri_from_spectator_nucleus_b_over_all_events = 0 ;
}


cell::~cell(){
}





