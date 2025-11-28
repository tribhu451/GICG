#include "mc_glau.h"

// #define ROTATED_SHIFT 

using std::cout;
using std::endl;

mc_glau::mc_glau(InputData *InData1)
{
  InData = InData1; 
  set_mc_glau_params();
  double gamma_k = InData->gamma_fluctuation_k ; 
    
  t1 = new random_gen();
  tr1 = new random_gen();
  f1 = new random_gen();
  f2 = new random_gen();
  w1 = new random_gen(gamma_k);
}

mc_glau::~mc_glau()
{
  delete t1;
  delete tr1;
  delete f1;
  delete f2;
  delete w1;
}


// This is the master function
void mc_glau::event()
{
  NPART = 1E5;
  NCOLL = 1E5;
  IMPACT_PARAM = 1E5;
  
  for(int j=0;j<=A;j++){XA[j]=0.0;YA[j]=0.0;ZA[j]=0.0;npart_tag_A[j]=0;}
  for(int j=0;j<=B;j++){XB[j]=0.0;YB[j]=0.0;ZB[j]=0.0;npart_tag_B[j]=0;}
  for(int j=0;j<500;j++){npart_x[j]=0.0;npart_y[j]=0.0;npart_w[j]=0.0;}
  for(int j=0;j<10000;j++){ncoll_x[j]=0.0;ncoll_y[j]=0.0;}
  
  for(int j=0;j<500;j++){npart_x_of_A[j]=0.0;npart_y_of_A[j]=0.0;npart_w_of_A[j]=0.0;}
  for(int j=0;j<500;j++){npart_x_of_B[j]=0.0;npart_y_of_B[j]=0.0;npart_w_of_B[j]=0.0;}

  shift_xavg_of_nucleons = 0. ;   
  shift_yavg_of_nucleons = 0. ;   
  
  //generate orientation angles of target & projectile ...
  double p_ori_theta = f2->random_sinx_0_to_pi();  
  double t_ori_theta = f2->random_sinx_0_to_pi();  
  double p_ori_phi = (2.0*M_PI)*(t1->rand_uniform());
  double t_ori_phi = (2.0*M_PI)*(t1->rand_uniform());
  //cout<<"[Info] (projectile orientation) p_theta: "<<p_ori_theta<<" p_phi: "<<p_ori_phi<<endl;
  //cout<<"[Info] (target orientation) t_theta: "<<t_ori_theta<<" t_phi: "<<t_ori_phi<<endl;
  
  //generate nucleus
  generate_nucleus(XA,YA,ZA,A,p_radius,p_dlt,p_beta2,p_beta4,p_ori_theta,p_ori_phi);
  generate_nucleus(XB,YB,ZB,B,t_radius,t_dlt,t_beta2,t_beta4,t_ori_theta,t_ori_phi);
  
  // generate impact parameter between bmin-bmax ...
  double b=f1->random_linear(bmin,bmax); 
  IMPACT_PARAM = b;                     
  //cout<<"[Info] b = "<<b<<" (fm)"<<endl;
  
  
  double zhi;
#ifdef ROTATED_SHIFT
  zhi=(2.0*M_PI)*(t1->Rndm());
#else
  zhi =0.0;
#endif
  
  //shifting of nucleus 
  shift_nucleus( XA, YA,  ZA, A, +b/2.0, zhi, XA, YA, ZA );
  shift_nucleus( XB, YB,  ZB, B, -b/2.0, zhi, XB, YB, ZB);
  
  
  // calculating npart & ncoll ...
  calculate_npart_ncoll(XA,YA,XB,YB,NPART,NCOLL,npart_x,npart_y,ncoll_x, ncoll_y);                    

}


void mc_glau::generate_nucleus(double* X1, double* Y1,double* Z1,int A,
			       double R, double dlt, double BETA2, double BETA4, double etaA, double psiA)
{    
  double X[500];double Y[500];double Z[500];
  
  double CMx=0.0;double CMy=0.0;double CMz=0.0;
  int count=0;
  
  do{
      double r=(15.0)*(tr1->rand_uniform());
      double Theta=(M_PI)*(tr1->rand_uniform());
      double Phi=((2.0)*M_PI)*(tr1->rand_uniform());
      double test=tr1->rand_uniform();
        
      double Y20=0.25*sqrt(5.0/M_PI)*
	(3*cos(Theta)*cos(Theta)-1.0);
      double Y40=(3.0/(16.0*sqrt(M_PI)))* 
	((35*pow(cos(Theta),4))-
	 (30*pow(cos(Theta),2))+3);
      double RAT= R*(1+(BETA2*Y20)+(BETA4*Y40));
      double rho=(1.0/60.0)*(r*r*(sin(Theta)))/(1.0+(exp((r-RAT)/dlt)));
      if( rho < 0 || rho > 1){
	 std::cout << "Problem in aaceptance rejection method" << std::endl ; 
         std::cout << " during generation of nuclues ..." << std::endl ; 
	 std::cout << "rho < 0 || rho > 1 :  rho = " << rho << " ..."<< std::endl ; 
	 exit(1);
      }      
      
      if(test < rho){
	  X[count]= (r*sin(Theta)*cos(Phi));
	  Y[count]=(r*sin(Theta)*sin(Phi));
	  Z[count]=(r*cos(Theta));
	  CMx=CMx+X[count]; CMy=CMy+Y[count] ;CMz=CMz+Z[count];    
	  count=count+1;
      }   
    }while(count<A);
  
  CMx=CMx/A;CMy=CMy/A;CMz=CMz/A;
  
  for(int j=0;j<A;j++){ X[j]=X[j]+(-CMx); Y[j]=Y[j]+(-CMy);Z[j]=Z[j]+(-CMz);}
  
  //etaA - nucleus orientaton angle (theta)
  //psiA - nucleus orientation angle (phi)
  
  for(int j=0;j<A;j++){
      X1[j]=(cos(psiA)*cos(etaA)*X[j])+(-sin(psiA)*Y[j])+(-cos(psiA)*sin(etaA)*Z[j]);
      Y1[j]=(sin(psiA)*cos(etaA)*X[j])+(cos(psiA)*Y[j])+(-sin(psiA)*sin(etaA)*Z[j]);
      Z1[j]=(sin(etaA)*X[j])+(cos(etaA)*Z[j]);
    }
}



// this function shifts the nucleus
void mc_glau::shift_nucleus(double* X1, double* Y1, double* Z1,int A, double b,
			    double zhi,double* X2, double* Y2, double* Z2 )
{
  for(int j=0;j<A;j++){
    X2[j]=X1[j]+((b)*cos(zhi));
    Y2[j]=Y1[j]+((b)*sin(zhi));
  }
}



// this function calculates N_{part} & N_{coll}
void mc_glau::calculate_npart_ncoll(double* vxA,double* vyA,double* vxB,double* vyB, int &Npart, 
 int &Ncoll, double* Npart_x, double* Npart_y, double* Ncoll_x, double* Ncoll_y){
  
  Ncoll                = 0 ;
  Npart                = 0 ;
  Nparticipants_from_A = 0 ; 
  Nparticipants_from_B = 0 ; 

  double weight;
  double occA[1000];
  double occB[1000];         //flag during calc of Npart
  //double Ncoll_x[2000]; double Ncoll_y[2000];  // x & y co-ordinate of binary collision sources
  //double Npart_x[1000]; double Npart_y[1000];  // x & y co-ordinate of participant sources
    
  for(int i=0;i<A;i++){
    occA[i]=0;
  }
  for(int i=0;i<B;i++){
    occB[i]=0;
  }
  
  for(int i=0; i<A; i++){
    for(int j=0; j<B; j++){  
	 double d=sqrt( pow((vxB[j]-vxA[i]),2) + 
			      pow ( (vyB[j]-vyA[i]),2));
	 double D=sqrt( ( sigma ) / (  M_PI ) ); 
	
	 if( d <= D){ 
	   Ncoll_x[Ncoll]=(vxA[i]+vxB[j])/2;
	   Ncoll_y[Ncoll]=(vyA[i]+vyB[j])/2;
	   Ncoll=Ncoll+1;
	  
	   if(occA[i]==0){
	    occA[i]=1;
         npart_tag_A[i] = 1 ; 
	    Npart_x[Npart]=vxA[i];
	    Npart_y[Npart]=vyA[i];
	    
	    weight = w1->rand_gamma();
	    npart_w[Npart]= weight;
	    
	    Npart=Npart+1;
	    npart_x_of_A[Nparticipants_from_A]=vxA[i];
	    npart_y_of_A[Nparticipants_from_A]=vyA[i]; 
	    npart_w_of_A[Nparticipants_from_A]=weight; 
         Nparticipants_from_A += 1 ;
	  } 
	  if(occB[j]==0){
	    occB[j]=1;
         npart_tag_B[j] = 1 ; 
	    Npart_x[Npart]=vxB[j];
	    Npart_y[Npart]=vyB[j];
	    
	    weight = w1->rand_gamma();
	    npart_w[Npart]=weight;
	    
	    Npart=Npart+1;
	    npart_x_of_B[Nparticipants_from_B]=vxB[j];
	    npart_y_of_B[Nparticipants_from_B]=vyB[j]; 
	    npart_w_of_B[Nparticipants_from_B]=weight; 
         Nparticipants_from_B += 1 ; 
	  } 
      }                                                           
     }                                                          
   }        
}



void mc_glau::get_nucleus_A(double *X1, double *Y1, double* Z1)
{
  for(int j=0;j<A;j++){
    X1[j]=XA[j]; 
    Y1[j]=YA[j];
    Z1[j]=ZA[j];
  }
}



void mc_glau::get_nucleus_B(double *X1, double *Y1, double* Z1)
{
  for(int j=0;j<B;j++){
    X1[j]=XB[j];
    Y1[j]=YB[j]; 
    Z1[j]=ZB[j]; 
  }
}


void mc_glau::get_npart_tag_in_nucleus_A(int *xx){
  for(int j=0;j<A;j++){
    xx[j] = npart_tag_A[j] ; 
  }
}


void mc_glau::get_npart_tag_in_nucleus_B(int *xx){
  for(int j=0;j<B;j++){
    xx[j] = npart_tag_B[j] ; 
  }
}




