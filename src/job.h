#pragma once
#include <fstream>
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <TRandom3.h>
#include <TF1.h>
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "input_data.h"
#include "mc_glau.h"
#include "cell.h"
#include "grid.h"
#include "idnni.h"
#include "mc_glau_smear.h"

#define PI  3.1415927

using std::cout;
using std::endl;
using std::to_string;


class job{

public :
 job(InputData* , grid* , mc_glau* );
 ~job();

 void event_by_event_mc_glauber_for_centrality_determination(int );
 void event_by_event_mc_glauber_and_eccentricity_calculation(int );

 void event_by_event_idnni_for_centrality_determination(int );
 void event_by_event_idnni_and_eccentricity_calculation(int );

 void rotate_by_second_order_participant_plane_and_then_gaussian_smearing_for_event_averaged_profile_using_idnni(int );
 void rotate_by_second_order_participant_plane_and_then_gaussian_smearing_for_event_averaged_profile_using_mc_glauber(int );

 void rotate_by_second_order_participant_plane_and_then_gaussian_smearing_for_ebe_mc_glauber_events(int event_no);

private :
 idnni* IDNNI ; 
 mc_glau* MC ; 
 InputData *iparams ;
 grid* arena ; 
 mc_glau_smear* MCGlbS ; 



};














