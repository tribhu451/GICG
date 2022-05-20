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

}

cell::~cell(){
}
