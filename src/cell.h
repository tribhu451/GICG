#pragma once

#include<fstream>
#include<string>
#include<sstream>
#include<iostream>
#include<fstream>
#include <TRandom3.h>
#include <TF1.h>
#include "TMath.h"
	
using std::cout;
using std::endl;

class cell{

 public :
  cell();
  ~cell();

 inline void set_position(int iix, int iiy) {
  ix = iix;
  iy = iiy;
 }

  // below functions are for each event calculationin IDNNI.
  inline void set_contri_from_nucleus_a(int xx){ contri_from_nucleus_a = xx ;}
  inline void set_contri_from_nucleus_b(int xx){ contri_from_nucleus_b = xx ;}

  inline double get_contri_from_nucleus_a(){return contri_from_nucleus_a;}
  inline double get_contri_from_nucleus_b(){return contri_from_nucleus_b;}
  inline double get_net_deposition_on_the_cell(){return (contri_from_nucleus_a+contri_from_nucleus_b);}

  inline void update_contri_from_nucleus_a(double xx){ contri_from_nucleus_a += xx;}
  inline void update_contri_from_nucleus_b(double xx){ contri_from_nucleus_b += xx;}
  

  // below functions are for gaussian smearing in 'each' event.  // NOTE : "each event."
  // event by event, one has to set the contribution inorder to start from fresh.  
  inline void set_contri_from_nucleus_a_after_gaussian_smearing(int xx){ contri_from_nucleus_a_after_gaussian_smearing = xx ;}
  inline void set_contri_from_nucleus_b_after_gaussian_smearing(int xx){ contri_from_nucleus_b_after_gaussian_smearing = xx ;}
  inline void set_contri_from_binary_collisions_after_gaussian_smearing(int xx){ contri_from_binary_collisions_after_gaussian_smearing = xx ;}

  inline double get_contri_from_nucleus_a_after_gaussian_smearing(){return contri_from_nucleus_a_after_gaussian_smearing;}
  inline double get_contri_from_nucleus_b_after_gaussian_smearing(){return contri_from_nucleus_b_after_gaussian_smearing;}
  inline double get_contri_from_binary_collisions_after_gaussian_smearing(){return contri_from_binary_collisions_after_gaussian_smearing;}
  inline double get_net_deposition_on_the_cell_after_gaussian_smearing(){
     return (contri_from_nucleus_a_after_gaussian_smearing+contri_from_nucleus_b_after_gaussian_smearing+contri_from_binary_collisions_after_gaussian_smearing);
  }

  inline void update_contri_from_nucleus_a_after_gaussian_smearing(double xx){ contri_from_nucleus_a_after_gaussian_smearing += xx;}
  inline void update_contri_from_nucleus_b_after_gaussian_smearing(double xx){ contri_from_nucleus_b_after_gaussian_smearing += xx;}
  inline void update_contri_from_binary_collisions_after_gaussian_smearing(double xx){ contri_from_binary_collisions_after_gaussian_smearing += xx;}


  // below functions are for event averaged profile.
  inline void update_contribution_from_nucleus_a_over_all_events(double xx){
    contri_from_nucleus_a_over_all_events += xx;
  }

  inline void update_contribution_from_nucleus_b_over_all_events(double xx){
    contri_from_nucleus_b_over_all_events += xx;
  }

  inline void update_contribution_from_binary_collisions_over_all_events(double xx){
    contri_from_binary_collisions_over_all_events += xx;
  }


  inline double get_contribution_from_nucleus_a_over_all_events(){
    return contri_from_nucleus_a_over_all_events ;
  }

  inline double get_contribution_from_nucleus_b_over_all_events(){
    return contri_from_nucleus_b_over_all_events  ;
  }

  inline double get_contribution_from_binary_collisions_over_all_events(){
    return contri_from_binary_collisions_over_all_events  ;
  }


  void reset_contributions_to_zero();


 private :

   int ix, iy;            // cell coordinate on the grid

   // contribution from nucleus A and Nucleus B.
   // Below to variables are to be used in IDNNI only.
   double  contri_from_nucleus_a;
   double  contri_from_nucleus_b;


   // below variables for gaussian smearing in each event.
   double  contri_from_nucleus_a_after_gaussian_smearing;
   double  contri_from_nucleus_b_after_gaussian_smearing;
   double  contri_from_binary_collisions_after_gaussian_smearing; // only required for 2 component MC Glauber. 


   // below variables are for event averaged profile.
   double  contri_from_nucleus_a_over_all_events ; 
   double  contri_from_nucleus_b_over_all_events ;
   double  contri_from_binary_collisions_over_all_events ; 





};
