#pragma once
#include <fstream>
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <TRandom3.h>
#include <TF1.h>
#include "TMath.h"
#include "input_data.h"
#include "mc_glau.h"
#include "cell.h"
#include <gsl/gsl_sf_erf.h>


using std::cout ; 
using std::endl ; 
using std::fstream ; 
using std::string ; 
using std::istringstream;


class rapidity_extension{
    public :
        rapidity_extension(InputData* );
        ~rapidity_extension();
        void write_rapidity_extended_tilted_profile_from_idnni_boost_invariant_deposition
                  (int flag_for_baryon, int flag_for_reading_external_file, std::string external_file_name );

    private :
        InputData *iparams;
        istringstream* iss;
        char   buff[400];
        double rapidity_asymmetry_ansatz_for_tilt_by_Bozek(double etas, double etam);
        int theta(double );
        double baryon_density_eta_envelop_profile(double eta, double eta0_nb, double sigma_eta_nb_plus, double sigma_eta_nb_minus ) ; 
        double baryon_density_eta_envelop_profile_function(double* x, double* p) ; 
        double integrate_baryon_density_eta_envelop_profile_over_eta() ; 

};



