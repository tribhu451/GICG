#pragma once
#include <fstream>
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include "input_data.h"
#include "mc_glau.h"
#include "cell.h"
#include "grid.h"
#include "mc_glau_smear.h"

#define PI  3.1415927

using std::cout;
using std::endl;
using std::to_string;
using std::istringstream;


class job{

public :
 job(InputData* , grid* , mc_glau* );
 ~job();


 void rotate_by_second_order_participant_plane_and_then_gaussian_smearing_for_event_averaged_profile_using_mc_glauber(int );
 void rotate_by_second_order_participant_plane_and_then_gaussian_smearing_for_ebe_mc_glauber_events(int event_no);

private :

 mc_glau* MC ; 
 InputData *iparams ;
 grid* arena ; 
 mc_glau_smear* MCGlbS ; 



};














