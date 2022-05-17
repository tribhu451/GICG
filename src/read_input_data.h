#include<iostream>
#include<fstream>
#include<sstream>
#include<string>
#include<chrono>
#include "input_data.h"

using std::cout;
using std::ofstream;
using std::endl;
using namespace std::chrono;
using std::string;
using std::cin;
using std::fstream;
using std::ios;
using std::istringstream;

class ReadInputPars{

public :

// This functions reads the input.dat file and sets the input parameters in the code 
void read_input_data_(InputData *input_parameter_list, string input_file_name)
{
  string a_; char a[50];
  
  istringstream* iss;
  char   buff[200];
  
  fstream File0;
  File0.open(input_file_name,ios::in);
  if(!File0){cout<<"No input file, exit..."; exit(1);}
  
  int number = 0;
  while(!File0.eof())
    {
      File0.getline(buff,200);
      if (!(*buff) || (*buff == '#')) {number ++; continue;}
      iss = new istringstream(buff);
      *iss >> a_ >> a ;
 
      if(a_ == "projectile" )    {input_parameter_list->projectile = a;}
      if(a_ == "target" )    {input_parameter_list->target = a;}
      if(a_ == "SNN" )    {input_parameter_list->SNN = atof(a);}
 
      
      if(a_ == "bmin" )  {input_parameter_list->bmin = atof(a);}       
      if(a_ == "bmax" )  {input_parameter_list->bmax = atof(a);}  
      
      if(a_ == "xhard" )  {input_parameter_list->xhard = atof(a);}       
      if(a_ == "npp" )  {input_parameter_list->npp = atof(a);} 

      if(a_ == "root_output_file_name" ) {input_parameter_list->root_output_file_name = a;}
      
      if(a_ == "xmax" )  {input_parameter_list->xmax = atof(a);}       
      if(a_ == "ymax" )  {input_parameter_list->ymax = atof(a);}       
      if(a_ == "nx" )    {input_parameter_list->nx = atof(a);}       
      if(a_ == "ny" )    {input_parameter_list->ny = atof(a);}       

      if(a_ == "gaussian_smearing_sigma" ) {input_parameter_list->gaussian_smearing_sigma = atof(a);}       

      if(a_ == "lower_mult_proxy_cut" ) {input_parameter_list->lower_mult_proxy_cut = atof(a);}       
      if(a_ == "upper_mult_proxy_cut" ) {input_parameter_list->upper_mult_proxy_cut = atof(a);} 


      // Rapidity extension
      if(a_ == "etamax" ) {input_parameter_list->etamax = atof(a);} 
      if(a_ == "neta" ) {input_parameter_list->neta = atof(a);} 
      if(a_ == "matter_eta_plateau" ) {input_parameter_list->matter_eta_plateau = atof(a);} 
      if(a_ == "matter_eta_fall" ) {input_parameter_list->matter_eta_fall = atof(a);} 
      if(a_ == "matter_tilt_param" ) {input_parameter_list->matter_tilt_param = atof(a);} 

      if(a_ == "baryon_rapidity_profile_eta_peak" ) {input_parameter_list->baryon_rapidity_profile_eta_peak = atof(a);} 
      if(a_ == "baryon_rapidity_profile_sigma_eta_plus" ) {input_parameter_list->baryon_rapidity_profile_sigma_eta_plus = atof(a);} 
      if(a_ == "baryon_rapidity_profile_sigma_eta_minus" ) {input_parameter_list->baryon_rapidity_profile_sigma_eta_minus = atof(a);} 

      delete iss;
      number++;
    } 
  
  File0.close();
  
}

 
};//class end

