#include <iostream>
#include <string>
#include <cmath>
#include <sstream>
#include "TFile.h"
#include "TTree.h"
#include "input_data.h"
#include "read_input_data.h"
#include "mc_glau.h"
#include "idnni.h"
#include "job.h"
#include "rap_extn.h"
#include "mc_glau_smear.h"

using std::cout;
using std::endl;
using std::to_string;
using std::istringstream;


int main(int argc, char **argv)
{
  // reading the input data from the file
  string input_file_name;char* event_no_s ;
  if(argc == 3){
    event_no_s = *(argv+1);
    input_file_name = *(argv+2);}
  else{
    cout<<"plz give 2 arguments\n" <<
          "1st argument - no of events you want to generate " << 
          " / the no. of configurations you want to put in order to produce avergaed profile.\n"
          "2nd argument - input filename"<<endl;
    exit(1);
   }
  ReadInputPars reader;
  InputData InData;
  reader.read_input_data_(&InData, input_file_name);
  int event_no = atof(event_no_s) ;


  cout << "\n\n" << endl ;
  cout<<"      *********************************************************************     "<<endl;
  cout<<"      *          GCIG(Glauber-type Initial Condition Generator)           *     "<<endl;
  cout<<"      *********************************************************************     "<<endl;
  cout << "\n\n" << endl ;

  cout<<"Total no of events : "<<event_no<<endl;  
  cout << InData.projectile << "+" << InData.target << " at " << InData.SNN << "GeV" << endl;

  grid* arena = new grid(&InData);
  mc_glau* MC = new mc_glau(&InData);
  job* JOB = new job(&InData,arena,MC); 
  rapidity_extension* RE = new rapidity_extension(&InData); 
 
  // ********** uncomment to perform jobs as per your requirement. ********** //
  //JOB->event_by_event_mc_glauber_for_centrality_determination(event_no);
  //JOB->event_by_event_mc_glauber_and_eccentricity_calculation(event_no);
  //JOB->event_by_event_idnni_for_centrality_determination(event_no);
  //JOB->event_by_event_idnni_and_eccentricity_calculation(event_no);
  //JOB->rotate_by_second_order_participant_plane_and_then_gaussian_smearing_for_event_averaged_profile_using_idnni(event_no);
  //JOB->rotate_by_second_order_participant_plane_and_then_gaussian_smearing_for_event_averaged_profile_using_mc_glauber(event_no); 
  //RE ->write_rapidity_extended_tilted_profile_from_idnni_boost_invariant_deposition(0,0,"xx"); // nb flag, external file flag, external file name ;
  //RE ->write_rapidity_extended_tilted_profile_from_mc_glauber_boost_invariant_deposition(1,0,"xx"); // nb flag, external file flag, external file name ;
  //RE ->write_rapidity_extended_energy_momentum_conserving_SARJ_profile_from_mc_glauber_boost_invariant_deposition(1,0,"xx"); 
  //RE ->write_rapidity_extended_Bjorken_flow_breaking_tilted_profile_from_mc_glauber_boost_invariant_deposition(1,0,"xx"); 

  JOB->rotate_by_second_order_participant_plane_and_then_gaussian_smearing_for_ebe_mc_glauber_events(event_no);

  for (int ii=0; ii<event_no; ii++){
    std::stringstream filename;
    filename.str("");
    filename << "output/mc_glauber_boost_invariant_event_averaged_profile_for_rapidity_extension_" << ii ;
    filename << ".dat";
    RE ->write_rapidity_extended_tilted_profile_from_mc_glauber_boost_invariant_deposition(1,1,filename.str().c_str(),ii); 
          // nb flag, external file flag, external file name ;
  }

  delete MC ;
  delete arena ;
  delete JOB ;
  delete RE ;    

  
  return 0;
}







