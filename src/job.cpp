#include "job.h"


job::job(InputData* iparams_, grid* arena_, mc_glau* MC_ ){
  iparams = iparams_ ; 
  arena = arena_ ; 
  MC = MC_ ;
  if(iparams->mode == 1){
    MCGlbS = new mc_glau_smear(iparams,arena,MC);
  }
  else{
    std::cout << "mode ==" << iparams->mode 
    << "  not available ... " << std::endl ;
    exit(-1);
  }
}


job::~job(){
}


void job::gaussian_smearing_for_ebe_mc_glauber_events(int event_no){
  std::ofstream outfile;
  std::stringstream output_filename;
  output_filename.str("");
  output_filename << "events_info.dat" ;
  outfile.open(output_filename.str().c_str(), std::ios::out);

  double b; int npart; double multiplicity_;
  double eccentricity_[6] = {0.}; 
  int event_count = 0 ;
  do{
     MC->event();
     b = MC->get_impactf();
     npart = MC->get_npart();
     if(npart<1){continue;}
     MCGlbS->smear_it(iparams->gaussian_smearing_sigma);
     multiplicity_ = MCGlbS->get_reduced_thickness_after_gaussian_smearing();
     MCGlbS->get_eccentricities_after_gaussian_smearing(eccentricity_);
     
     //print
     outfile << event_count << "  " << b << "  " << npart << "   " << multiplicity_ << "  " 
     << eccentricity_[2] << "  " << eccentricity_[3] << "  " << eccentricity_[4] 
     << "  " << eccentricity_[5] << std::endl ;
      
     // write to files
     //MCGlbS->write_event_averaged_profile_to_file_after_gaussian_smearing(1,1,event_count);
     //MCGlbS->write_event_averaged_profile_to_file_after_gaussian_smearing(1,0,event_count);
     
     // reset everything and cleanup for next event
     MCGlbS->reset_contribution_to_zero_on_the_cells();
     event_count++;
  }while( event_count < event_no ) ;

  outfile.close();
}






