#include "job.h"


job::job(InputData* iparams_, grid* arena_, mc_glau* MC_ ){
  iparams = iparams_ ; 
  arena = arena_ ; 
  MC = MC_ ;
  if(iparams->mode == 1){
    MCGlbS = new mc_glau_smear(iparams,arena,MC);
  }
  else{
   std::cout << "mode ==" << iparams->mode << "  X " << std::endl ; exit(-1);
  }
}


job::~job(){
}



void job::rotate_by_second_order_participant_plane_and_then_gaussian_smearing_for_event_averaged_profile_using_mc_glauber(int event_no){

  double b;
  int event_count = 0 ;
  do 
    {
      // It is must in all cases to generate nucleons.
      MC->event(0);
      b = MC->get_impactf();
      if( MC->get_two_component_galuber_multiplicity_proxy() > iparams->lower_mult_proxy_cut && 
             MC->get_two_component_galuber_multiplicity_proxy() <= iparams->upper_mult_proxy_cut ){

           if( (event_count%50) == 0 ){
             std::cout << "event no = " << event_count << ",  b = " << b << " (fm),  Npart = " << MC->get_npart() << std::endl ; 
           }
           MCGlbS->smear_it(iparams->gaussian_smearing_sigma);
           MCGlbS->update_contribution_on_cells_over_all_events_with_gaussian_smearing();
           event_count++;
      }
    }
  while( event_count < event_no ) ;

  // [About below 2 lines] The second argument is a flag for making a boost invariant profile for music input
  // The third argument is for the event index, which is redundant if second argument is 1.
  MCGlbS->write_event_averaged_profile_to_file_after_gaussian_smearing(event_no,0,0); 
  MCGlbS->write_event_averaged_profile_to_file_after_gaussian_smearing(event_no,1,0);

}


void job::rotate_by_second_order_participant_plane_and_then_gaussian_smearing_for_ebe_mc_glauber_events(int event_no){

  std::ofstream outfile;
  std::stringstream output_filename;

  double b;
  int event_count = 0 ;
  do 
    {
      // It is must in all cases to generate nucleons.
      // Here we are also calculating eccentricities and
      // participant angles //
      MC->event(1);
      b = MC->get_impactf();
      if( MC->get_two_component_galuber_multiplicity_proxy() > iparams->lower_mult_proxy_cut && 
             MC->get_two_component_galuber_multiplicity_proxy() <= iparams->upper_mult_proxy_cut ){
           std::cout << "event no = " << event_count << ",  b = " << b << " (fm),  Npart = " << MC->get_npart() << std::endl ; 
           MCGlbS->smear_it(iparams->gaussian_smearing_sigma);
           MCGlbS->update_contribution_on_cells_over_all_events_with_gaussian_smearing();
           MCGlbS->write_event_averaged_profile_to_file_after_gaussian_smearing(1,1,event_count);
           MCGlbS->write_event_averaged_profile_to_file_after_gaussian_smearing(1,0,event_count);
           MCGlbS->reset_contribution_from_all_events_to_zero_on_the_cells();
           event_count++;
      }
    }
  while( event_count < event_no ) ;


}






