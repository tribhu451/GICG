#include <iostream>
#include <string>
#include <cmath>
#include <sstream>
#include "input_data.h"
#include "read_input_data.h"
#include "mc_glau.h"
#include "job.h"
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

  JOB->gaussian_smearing_for_ebe_mc_glauber_events(event_no); 
 
  delete MC ;
  delete arena ;
  delete JOB ;
  
  return 0;
}







