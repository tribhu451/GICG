// *************************** //
// *  25 feb. 2021           * //
// *  tribhu.451@gmail.com   * //
// *  version 2.0            * //
// *************************** //


//  This class calculates npart, ncoll and eccentricity
//  in mc glauber model

#pragma once

#include<fstream>
#include<string>
#include<sstream>
#include<iostream>
#include<fstream>
#include <TRandom3.h>
#include <TF1.h>
#include "TMath.h"
#include "input_data.h"
	
using std::cout;
using std::endl;

class mc_glau
{
  
 public:

  mc_glau(InputData *InData1);
  ~mc_glau();
  void event(int flag_for_centrality_calculation );
  inline int get_npart(){return NPART;}
  inline int get_no_of_participants_in_nucleus_a(){return Nparticipants_from_A;}
  inline int get_no_of_participants_in_nucleus_b(){return Nparticipants_from_B;}
  inline int get_no_of_spectators_in_nucleus_a(){return spectators_from_A;}
  inline int get_no_of_spectators_in_nucleus_b(){return spectators_from_B;}
  inline int get_ncoll(){return NCOLL;}
  inline double get_two_component_galuber_multiplicity_proxy(){
    return  InData->npp * 0.5 * get_npart()* ( 1.0 - InData->xhard ) +  InData->npp * InData->xhard * get_ncoll() ; 
  }
  inline double get_impactf(){return IMPACT_PARAM;}
  void calculate_eccentricity(int Norder, int aN_part,int aN_coll,double *Npart_x,
                     double *Npart_y,double *Ncoll_x,double *Ncoll_y, double* eps, double* psi);

  void get_nucleus_A(double *X1, double *Y1, double* Z1, int* );
  void get_nucleus_B(double *X2, double *Y2, double* Z2, int* );

  inline int get_mass_number_of_nucleus_A(){return A;}
  inline int get_mass_number_of_nucleus_B(){return B;}
  inline int get_atomic_number_of_nucleus_A(){return AZ;}
  inline int get_atomic_number_of_nucleus_B(){return BZ;}


  inline double get_radius_of_nucleus_A(){return p_radius;}
  inline double get_radius_of_nucleus_B(){return t_radius;}

  inline double get_inelastic_cross_section(){return sigma; }

  inline double eccen1(){return eccentricity[1];}
  inline double phi1(){return PhiN[1];}
  inline double eccen2(){return eccentricity[2];}
  inline double phi2(){return PhiN[2];}
  inline double eccen3(){return eccentricity[3];}
  inline double phi3(){return PhiN[3];}
  inline double eccen4(){return eccentricity[4];}
  inline double phi4(){return PhiN[4];}
  inline double eccen5(){return eccentricity[5];}
  inline double phi5(){return PhiN[5];}
  inline double eccen6(){return eccentricity[6];}
  inline double phi6(){return PhiN[6];}

  inline void get_npart_source_positions(double *xx, double* yy){
    for(int ii=0; ii<500; ii++){
      xx[ii] = npart_x[ii] ; 
      yy[ii] = npart_y[ii] ; 
    }
  }

  inline void get_ncoll_source_positions(double *xx, double* yy){
    for(int ii=0; ii<10000; ii++){
      xx[ii] = ncoll_x[ii] ; 
      yy[ii] = ncoll_y[ii] ; 
    }
  }


  inline void get_npart_source_positions_of_nucleus_a(double *xx, double* yy){
    for(int ii=0; ii<500; ii++){
      xx[ii] = npart_x_of_A[ii] ; 
      yy[ii] = npart_y_of_A[ii] ; 
    }
  }


  inline void get_npart_source_positions_of_nucleus_b(double *xx, double* yy){
    for(int ii=0; ii<500; ii++){
      xx[ii] = npart_x_of_B[ii] ; 
      yy[ii] = npart_y_of_B[ii] ; 
    }
  }


  inline void get_spectator_positions_and_its_proton_flags_of_nucleus_a(double *xx, double* yy, int* pflag){
    for(int ii=0; ii<500; ii++){
      xx[ii] = spec_x_of_A[ii] ; 
      yy[ii] = spec_y_of_A[ii] ; 
      pflag[ii] = spec_proton_flags_of_A[ii] ; 
    }
  }


  inline void get_spectator_positions_and_its_proton_flags_of_nucleus_b(double *xx, double* yy, int* pflag){
    for(int ii=0; ii<500; ii++){
      xx[ii] = spec_x_of_B[ii] ; 
      yy[ii] = spec_y_of_B[ii] ; 
      pflag[ii] = spec_proton_flags_of_B[ii] ; 
    }
  }

 
 private:
  InputData *InData;
  
  int A; // mass no. of projrctile nucleus
  int B; // mass no. of target nucleus
  int AZ; // atomic no. of projrctile nucleus
  int BZ; // atomic no. of target nucleus
  double sigma; // energy-> cross-section
  
// projectie wood-Saxon parameters
  double p_radius;        
  double p_dlt;        
  double p_beta2;
  double p_beta4;

// target wood-Saxon parameters
  double t_radius;        
  double t_dlt;        
  double t_beta2;
  double t_beta4;

// npart, ncoll and impact parameter
  int NPART;
  int NCOLL;
  double IMPACT_PARAM;

  // individual contribution
  int Nparticipants_from_A ; 
  int Nparticipants_from_B ; 
  int spectators_from_A ; 
  int spectators_from_B ; 


// two-component energy deposition
  double npp;
  double X_hard;

 //impact parameter range
  double bmin,bmax;

// eccentricity and partcipant plane angle
  double eccentricity[10];
  double PhiN[10];

  TRandom3* t1;
  TRandom3* tr1;
  TRandom3* tr2;
  TF1* f1;
  TF1* f2;

  // XA, YA, ZA -> (x,y,z) coordinate of nucleons in nucleus A.
  // XB, YB, ZB -> (x,y,z) coordinate of nucleons in nucleus B.
  // PROT_FLAG_X -> indicates wheather the nucleon is a proton or not. 
  double XA[300];double YA[300];double ZA[300]; int PROT_FLAG_A[300];
  double XB[300];double YB[300];double ZB[300]; int PROT_FLAG_B[300];  
  double npart_x[500],npart_y[500];
  double ncoll_x[10000],ncoll_y[10000];

  
  // x and y coordinate of nparts coming from nucleus A and B.
  double npart_x_of_A[500],npart_y_of_A[500];
  double npart_x_of_B[500],npart_y_of_B[500];

  // x and y cordinate of spectators and the flag to indicate 
  // whether the spectator is a proton or not.
  // "A" stands for projectile whereas "B" stands for target nucleus.
  double spec_x_of_A[500],spec_y_of_A[500],spec_proton_flags_of_A[500];
  double spec_x_of_B[500],spec_y_of_B[500],spec_proton_flags_of_B[500];



  void generate_nucleus(double* X1, double* Y1,double* Z1,int* PROT_FLAG___, int A, int ATOMIC_NO,
			double R, double dlt, double BETA2, double BETA4, double etaA, double psiA);
  
  void calculate_npart_ncoll(double* vxA,double* vyA, int* , double* vxB,double* vyB, int* , int &Npart, 
			     int &Ncoll, double* Npart_x, double* Npart_y, double* Ncoll_x, double* Ncoll_y);
  
  void shift_nucleus(double* X1, double* Y1, double* Z1,int A, double b,
                     double zhi,double* X2, double* Y2, double* Z2 );





void set_mc_glau_params()
{
  npp = InData->npp;
  X_hard = InData->xhard;

  // projectile nucleus
  if(InData->projectile == "Au") {A=197; AZ=79; p_radius = 6.42; p_dlt=0.41; p_beta2 = -0.13; p_beta4 =0.0;} //arXiv: 1409.8375 [Table. 1]
  else if(InData->projectile == "Au2") {A=197; AZ=79; p_radius = 6.37; p_dlt=0.53; p_beta2 = 0.0; p_beta4 =0.0;} //arXiv: 1409.8375 [Table. 1]
  else if(InData->projectile == "Pb") {A=208; AZ=82; p_radius = 6.66; p_dlt=0.45; p_beta2 = 0.0; p_beta4 =0.0;}
  else if(InData->projectile == "p") {A=1; AZ=1; p_radius = 0.1; p_dlt=0.1; p_beta2 = 0.0; p_beta4 =0.0;} // woods-saxon params can be ignored-
  //because we will sample only one nucleon.
  else if(InData->projectile == "U") {A=238; AZ=92; p_radius = 6.86; p_dlt=0.42; p_beta2 = 0.265; p_beta4 =0.093;}
  else if(InData->projectile == "U2") {A=238; AZ=92; p_radius = 6.86; p_dlt=0.42; p_beta2 = 0.0; p_beta4 =0.0;}
  else {cout<<"projectile not recognized, it's : "<<InData->projectile<<endl; exit(1);}
  
  // target nucleus
  if(InData->target == "Au") {B=197; BZ=79; t_radius = 6.42; t_dlt=0.41; t_beta2 = -0.13; t_beta4 =0.0;}
  else if(InData->target == "Au2") {B=197; BZ=79; t_radius = 6.37; t_dlt=0.53; t_beta2 = 0.00; t_beta4 =0.0;} 
  else if(InData->target == "Pb") {B=208; BZ=82, t_radius = 6.66; t_dlt=0.45; t_beta2 = 0.0; t_beta4 =0.0;}
  else if(InData->target == "p") {B=1; BZ=1; t_radius = 0.1; t_dlt=0.1; t_beta2 = 0.0; t_beta4 =0.0;}
  else if(InData->target == "U"){B=238;  BZ=92; t_radius = 6.86; t_dlt=0.42; t_beta2 = 0.265; t_beta4 =0.093;}
  else if(InData->target == "U2"){B=238; BZ=92; t_radius = 6.86; t_dlt=0.42; t_beta2 = 0.0; t_beta4 =0.0;}
  else {cout<<"target not recognized, it's : "<<InData->target<<endl; exit(1);}
  

  // collision energy and inelastic cross-section
  if (InData->SNN > 5.0 ){
    //  https://arxiv.org/pdf/1901.04484.pdf ( Equation-4 )  //
    double ecm = InData->SNN ; 
    double sigma_NN_inel  = 40.32 * pow( ecm + 53.08, 0.104) - 30.15 - 8.75 / ecm ; 
    sigma = (sigma_NN_inel*0.1) ;
  }
  //  https://arxiv.org/pdf/1712.05815.pdf (Table 4.2 page-106)  //
  //  https://arxiv.org/pdf/nucl-ex/0701025.pdf (caption of figure-2)  //
  if (InData->SNN == 200.0 )  {  sigma =4.20;   }
  if (InData->SNN == 130.0 )  {  sigma =4.00;   }
  if (InData->SNN == 62.4  )  {  sigma =3.56;   }
  if (InData->SNN == 19.6  )  {  sigma =3.23;   }
  if (InData->SNN == 2760.0)  {  sigma =6.40;   }
  if (InData->SNN == 5500.0)  {  sigma =7.20;   }
  if (InData->SNN < 5.0 ){ 
     cout<<"Glauber inelastic crossection not available ..."<<endl; 
     exit(1);
  }
  
  // impact parameter range
  bmin = InData->bmin; bmax = InData->bmax;
  
  if(InData->projectile != "p")
    {
      cout<<"[Info] projectile mass number : "<<A<<endl;
      cout<<"[Info] wood-saxon | projectile Radius : "<<p_radius<<"   a: "<<p_dlt<<endl;
      cout<<"[Info] deformation parameter of projectile | Beta-2 : "<<p_beta2<<" and Beta-4 : "<<p_beta4<<"\n"<<endl; 
    }
  
  if(InData->target != "p")
    {
      cout<<"[Info] target mass number : "<<B<<endl;
      cout<<"[Info] wood-saxon | target Radius : "<<t_radius<<"   a: "<<t_dlt<<endl;
      cout<<"[Info] deformation parameter of target | Beta-2 : "<<t_beta2<<" and Beta-4 : "<<t_beta4<<"\n"<<endl;
    }
  
}



};
















