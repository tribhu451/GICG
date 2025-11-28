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


void job::gaussian_smearing_for_ebe_mc_glauber_events(int event_no){

  std::ofstream outfile;
  std::stringstream output_filename;

  double b; int npart;
  int event_count = 0 ;
  do 
    {
      // It is must in all cases to generate nucleons.
      // Here we are also calculating eccentricities and
      // participant angles //
      MC->event();
      b = MC->get_impactf();
      npart = MC->get_npart();
      std::cout << event_count << "  " << b << "  " << npart << "    0 " << std::endl ; 
      //MCGlbS->smear_it(iparams->gaussian_smearing_sigma);
      //MCGlbS->update_contribution_on_cells_over_all_events_with_gaussian_smearing();
      //MCGlbS->write_event_averaged_profile_to_file_after_gaussian_smearing(1,1,event_count);
      //MCGlbS->write_event_averaged_profile_to_file_after_gaussian_smearing(1,0,event_count);
      //MCGlbS->reset_contribution_from_all_events_to_zero_on_the_cells();
      event_count++;
    }
  while( event_count < event_no ) ;


}






