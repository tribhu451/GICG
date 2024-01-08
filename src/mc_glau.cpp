#include "mc_glau.h"

// #define ROTATED_SHIFT 

using std::cout;
using std::endl;

mc_glau::mc_glau(InputData *InData1)
{
  InData = InData1; 
  set_mc_glau_params();
    
  // Random number generators to be used in function event().
  t1=new TRandom3();
  t1->SetSeed(0);
  long kss=t1->GetSeed();
  gRandom->SetSeed(kss);
  f1= new TF1("f1","x",0.0,25.0);
  f2= new TF1("f2","sin(x)",0.0,TMath::Pi());
  
  // Random number generators to be used in function generate_nucleus().
  tr1 = new TRandom3();
  tr1->SetSeed(0);
  tr2 = new TRandom3();
  tr2->SetSeed(0);

}

mc_glau::~mc_glau()
{
  delete t1;
  delete tr1;
  delete f1;
  delete f2;
}


// This is the master function
void mc_glau::event(int flag_for_eccentricity_calculation)
{
  NPART = 1E5;
  NCOLL = 1E5;
  IMPACT_PARAM = 1E5;
  
  for(int j=0;j<=A;j++){XA[j]=1000.;YA[j]=1000.;ZA[j]=1000.;PROT_FLAG_A[j]=0;}
  for(int j=0;j<=B;j++){XB[j]=1000.;YB[j]=1000.;ZB[j]=1000.;PROT_FLAG_B[j]=0;}
  for(int j=0;j<500;j++){npart_x[j]=1000.;npart_y[j]=1000.;}
  for(int j=0;j<10000;j++){ncoll_x[j]=1000.;ncoll_y[j]=1000.;}
  
  for(int j=0;j<500;j++){npart_x_of_A[j]=1000.;npart_y_of_A[j]=1000.;}
  for(int j=0;j<500;j++){npart_x_of_B[j]=1000.;npart_y_of_B[j]=1000.;}

  for(int j=0;j<500;j++){spec_x_of_A[j]=1000.;spec_y_of_A[j]=1000.;spec_proton_flags_of_A[j]=0;}
  for(int j=0;j<500;j++){spec_x_of_B[j]=1000.;spec_y_of_B[j]=1000.;spec_proton_flags_of_B[j]=0;}

  
  //generate orientation angles of target & projectile ...
  double p_ori_theta = f2->GetRandom(0.0,TMath::Pi());  
  double t_ori_theta = f2->GetRandom(0.0,TMath::Pi());  
  double p_ori_phi = (2.0*TMath::Pi())*(t1->Rndm());
  double t_ori_phi = (2.0*TMath::Pi())*(t1->Rndm());
  //cout<<"[Info] (projectile orientation) p_theta: "<<p_ori_theta<<" p_phi: "<<p_ori_phi<<endl;
  //cout<<"[Info] (target orientation) t_theta: "<<t_ori_theta<<" t_phi: "<<t_ori_phi<<endl;
  
  //generate nucleus
  generate_nucleus(XA,YA,ZA,PROT_FLAG_A,A,AZ,p_radius,p_dlt,p_beta2,p_beta4,p_ori_theta,p_ori_phi);
  generate_nucleus(XB,YB,ZB,PROT_FLAG_B,B,BZ,t_radius,t_dlt,t_beta2,t_beta4,t_ori_theta,t_ori_phi);
  
  // generate impact parameter between bmin-bmax ...
  double b=f1->GetRandom(bmin,bmax); 
  IMPACT_PARAM = b;                     
  //cout<<"[Info] b = "<<b<<" (fm)"<<endl;
  
  
  double zhi;
#ifdef ROTATED_SHIFT
  zhi=(2.0*TMath::Pi())*(t1->Rndm());
#else
  zhi =0.0;
#endif
  
  //shifting of nucleus 
  shift_nucleus( XA, YA,  ZA, A, +b/2.0, zhi, XA, YA, ZA );
  shift_nucleus( XB, YB,  ZB, B, -b/2.0, zhi, XB, YB, ZB);
  
  
  // calculating npart & ncoll ...
  calculate_npart_ncoll(XA,YA,PROT_FLAG_A,XB,YB,PROT_FLAG_B,NPART,NCOLL,npart_x,npart_y,ncoll_x, ncoll_y);                    
  //cout<<"[Info] No. of participants : "<<NPART<<endl;
  //cout<<"[Info] No. of binary collisions : "<<NCOLL<<endl;
  
  
  // calculating eccentricity ...
  if (flag_for_eccentricity_calculation > 0 ){
    int Norder = 7 ;
    double epspp[Norder];
    double phipp[Norder];
    for(int iorder = 0 ; iorder < Norder ; iorder++ ){
      epspp[iorder] = 0. ; 
      phipp[iorder] = 0. ; 
    }
      calculate_eccentricity(Norder,NPART,NCOLL,npart_x,npart_y,ncoll_x,ncoll_y,epspp,phipp);
      for(int iorder = 0 ; iorder < Norder ; iorder++ ){
        eccentricity[iorder] = epspp[iorder] ; 
        PhiN[iorder]         = phipp[iorder] ; 
      }

  }
  
}


void mc_glau::generate_nucleus(double* X1, double* Y1,double* Z1,int* PROT_FLAG_, int A, int ATOMIC_NO,
			       double R, double dlt, double BETA2, double BETA4, double etaA, double psiA)
{    
  double X[500];double Y[500];double Z[500];
  
  double CMx=0.0;double CMy=0.0;double CMz=0.0;
  int count=0;
  
  do
    {
      double r=(11.0)*(tr1->Rndm());
      double Theta=(TMath::Pi())*(tr1->Rndm());
      double Phi=((2.0)*TMath::Pi())*(tr1->Rndm());
      double test=tr1->Rndm();
      
      
      double Y20=0.25*TMath::Sqrt(5.0/TMath::Pi())*
	(3*TMath::Cos(Theta)*TMath::Cos(Theta)-1.0);
      double Y40=(3.0/(16.0*TMath::Sqrt(TMath::Pi())))* 
	((35*TMath::Power(TMath::Cos(Theta),4))-
	 (30*TMath::Power(TMath::Cos(Theta),2))+3);
      double RAT= R*(1+(BETA2*Y20)+(BETA4*Y40));
      double rho=(1.0/60.0)*(r*r*(TMath::Sin(Theta)))/(1.0+(TMath::Exp((r-RAT)/dlt)));

      if(rho>1.0){
        std::cout << "[generate nucleus func in MC] : nuclear" <<
            " density (rho) > 1 => Not normalized ..." << std::endl;
        exit(1);
      }
      
      if(test < rho )
	{      
	  X[count]=(r*TMath::Sin(Theta)*TMath::Cos(Phi));
	  Y[count]=(r*TMath::Sin(Theta)*TMath::Sin(Phi));
	  Z[count]=(r*TMath::Cos(Theta));
	  CMx=CMx+X[count]; CMy=CMy+Y[count] ;CMz=CMz+Z[count];    
	  count=count+1;
	}   
    }   
  while(count<A);
  
  CMx=CMx/A;CMy=CMy/A;CMz=CMz/A;
  
  for(int j=0;j<A;j++){ X[j]=X[j]+(-CMx); Y[j]=Y[j]+(-CMy);Z[j]=Z[j]+(-CMz);}
  
  //etaA - nucleus orientaton angle (theta)
  //psiA - nucleus orientation angle (phi)
  
  for(int j=0;j<A;j++)
    {
      X1[j]=(TMath::Cos(psiA)*TMath::Cos(etaA)*X[j])+(-TMath::Sin(psiA)*Y[j])+(-TMath::Cos(psiA)*TMath::Sin(etaA)*Z[j]);
      Y1[j]=(TMath::Sin(psiA)*TMath::Cos(etaA)*X[j])+(TMath::Cos(psiA)*Y[j])+(-TMath::Sin(psiA)*TMath::Sin(etaA)*Z[j]);
      Z1[j]=(TMath::Sin(etaA)*X[j])+(TMath::Cos(etaA)*Z[j]);
    }


  // set proton flag to each nucleon //
  /*
  int prot_counter = 0 ; 
  double random_ ; 
  do{
   for(int ii=0;ii<A;ii++){
     if(prot_counter<ATOMIC_NO){
       random_ = tr2->Rndm();
       if(random_<0.5 && PROT_FLAG[ii]==0){
         PROT_FLAG[ii] =1 ; 
         prot_counter +=1;
       }
     }
   }
  }while(prot_counter<ATOMIC_NO);

  if(prot_counter==ATOMIC_NO){
    // 
  }
  else{
    std::cout << "Desired amount of proton flag are" 
          << "not assigned to nucleons ..." << std::endl ;  
    exit(1);
  }
  */

  // Sandeep suggested to assign first Z(atomic no.) nucleons as protons.
  // 23rd oct. 2023 
  for(int ii=0;ii<ATOMIC_NO;ii++){
    PROT_FLAG_[ii] = 1 ; 
  }


}



// this function shifts the nucleus
void mc_glau::shift_nucleus(double* X1, double* Y1, double* Z1,int A, double b,
			    double zhi,double* X2, double* Y2, double* Z2 )
{
  for(int j=0;j<A;j++){
    X2[j]=X1[j]+((b)*TMath::Cos(zhi));
    Y2[j]=Y1[j]+((b)*TMath::Sin(zhi));
  }
}



// This function calculates N_{part} & N_{coll}.
// [Important] This function is only correct if we want to deposit energy according to two-component Glauber model.
// In this function, the participant and binary collsion sources are shifted by a vector (xaverage, yaverage) in 
// the transverse plane in such a way that the centre of mass of the energy profile remains in (x=0,y=0). 
// The spectator sources are also shifted by the same amount.
void mc_glau::calculate_npart_ncoll(double* vxA,double* vyA,int* pflagA, double* vxB,double* vyB, int* pflagB, int &Npart, 
				    int &Ncoll, double* Npart_x, double* Npart_y, double* Ncoll_x, double* Ncoll_y)
{
  
  Ncoll                = 0 ;
  Npart                = 0 ;
  Nparticipants_from_A = 0 ; 
  Nparticipants_from_B = 0 ; 


  int occA[1000];
  int occB[1000];         //flag during calc of Npart
    
  for(int i=0;i<A;i++){
    occA[i]=0;
  }
  for(int i=0;i<B;i++){
    occB[i]=0;
  }
  
  for (int i=0; i<A; i++){
      for (int j=0; j<B; j++){  
	double d=TMath::Sqrt( TMath::Power((vxB[j]-vxA[i]),2) + 
			      TMath::Power ( (vyB[j]-vyA[i]),2));
	double D=TMath::Sqrt( ( sigma ) / (  TMath::Pi() ) ); 
	
	if( d <= D){ 
	  Ncoll_x[Ncoll]=(vxA[i]+vxB[j])/2;
	  Ncoll_y[Ncoll]=(vyA[i]+vyB[j])/2;
	  Ncoll=Ncoll+1;
	  
	  if(occA[i]==0){ 
	    occA[i]=1;
	    Npart_x[Npart]=vxA[i];
	    Npart_y[Npart]=vyA[i];
	    Npart=Npart+1;
	    npart_x_of_A[Nparticipants_from_A]=vxA[i];
	    npart_y_of_A[Nparticipants_from_A]=vyA[i]; 
            Nparticipants_from_A += 1 ;
	  } 
	  if(occB[j]==0){
	    occB[j]=1;
	    Npart_x[Npart]=vxB[j];
	    Npart_y[Npart]=vyB[j];
	    Npart=Npart+1;
	    npart_x_of_B[Nparticipants_from_B]=vxB[j];
	    npart_y_of_B[Nparticipants_from_B]=vyB[j]; 
            Nparticipants_from_B += 1 ; 
	  }
	  
	}                                                           
      }                                                          
  } 

  // determine the spectator positions //
  spectators_from_A = 0 ; 
  spectators_from_B = 0 ;
  // spectators from A 
  for(int ii=0; ii<A; ii++){
    if(occA[ii]>0){
       continue ;
    }
    spec_x_of_A[spectators_from_A] = vxA[ii] ;   
    spec_y_of_A[spectators_from_A] = vyA[ii] ;   
    spec_proton_flags_of_A[spectators_from_A] = pflagA[ii] ;   
    spectators_from_A += 1 ; 
  }
  // spectators from B
  for(int jj=0; jj<B; jj++){
    if(occB[jj]>0){
       continue ;
    }
    spec_x_of_B[spectators_from_B] = vxB[jj] ;   
    spec_y_of_B[spectators_from_B] = vyB[jj] ;  
    spec_proton_flags_of_B[spectators_from_B] = pflagB[jj] ;    
    spectators_from_B += 1 ; 
  }

                                                        
  // [Info]  shifting the energy distributions center to (0,0) //   
  double xref1=0.0;
  double yref1=0.0;
  double wref1=0.0;
  for(int k=0;k<Npart;k++){ 
    xref1=xref1+(Npart_x[k]*(0.5*npp*(1-X_hard) ));
    yref1=yref1+(Npart_y[k]*(0.5*npp*(1-X_hard)));
    wref1=wref1+(0.5*npp*(1-X_hard));
  }
  
  double xref2=0.0;
  double yref2=0.0;
  double wref2=0.0;
  for(int k=0;k<Ncoll;k++){ 
    xref2=xref2+(Ncoll_x[k]*(npp*(X_hard)));
    yref2=yref2+(Ncoll_y[k]*(npp*(X_hard)));
    wref2=wref2+(npp*(X_hard));
  }
  
  double xAverage=((xref1+xref2)/(wref1+wref2));
  double yAverage=((yref1+yref2)/(wref1+wref2));
  
  // cout<<xAverage<<"  "<<yAverage<<"\n";
  
  // shift all the Npart, Ncoll and spectator sources by (xAverage,yAverage) //
  for(int k=0;k<Npart;k++){ 
    Npart_x[k] = Npart_x[k]-xAverage;
    Npart_y[k] = Npart_y[k]-yAverage;
  }
  for(int k=0;k<Ncoll;k++){ 
    Ncoll_x[k] = Ncoll_x[k]-xAverage;
    Ncoll_y[k] = Ncoll_y[k]-yAverage;
  }
  for(int k=0;k<Nparticipants_from_A;k++){ 
    npart_x_of_A[k] = npart_x_of_A[k]-xAverage;
    npart_y_of_A[k] = npart_y_of_A[k]-yAverage;
  }
  for(int k=0;k<Nparticipants_from_B;k++){ 
    npart_x_of_B[k] = npart_x_of_B[k]-xAverage;
    npart_y_of_B[k] = npart_y_of_B[k]-yAverage;
  }
  for(int ii=0;ii<spectators_from_A;ii++){ 
    spec_x_of_A[ii] = spec_x_of_A[ii]-xAverage;
    spec_y_of_A[ii] = spec_y_of_A[ii]-yAverage;
  }
  for(int jj=0;jj<spectators_from_B;jj++){ 
    spec_x_of_B[jj] = spec_x_of_B[jj]-xAverage;
    spec_y_of_B[jj] = spec_y_of_B[jj]-yAverage;
  }

}


// This function calculates eccentricity and participant plane angle
void mc_glau::calculate_eccentricity(int Norder, int aN_part,int aN_coll,double *Npart_x,
				     double *Npart_y,double *Ncoll_x,double *Ncoll_y,double* eps, double* psi)
{
  
  for(int i=0; i<Norder; i++){
    eps[i] = 0.0;
    psi[i] = 0.0;
  }
  
  int Total_Nch=aN_part+aN_coll;
  double Nch_r[Total_Nch];
  double Nch_phi[Total_Nch];
  double Nch_value[Total_Nch];
  
  double Npart_r[aN_part];
  double Npart_phi[aN_part];
  
  double Ncoll_r[aN_coll];
  double Ncoll_phi[aN_coll];
  
  double Ra,Theta;
  
  for(int k=0;k<aN_part;k++){ 
    Ra= TMath::Sqrt( TMath::Power(Npart_x[k],2)+TMath::Power( Npart_y[k],2) );
    Theta =(TMath::ATan2(Npart_y[k],Npart_x[k]));
    Npart_r[k]=Ra;
    Npart_phi[k]=Theta;
  }
  for(int k=0;k<aN_coll;k++){      
    Ra=TMath::Sqrt( TMath::Power(Ncoll_x[k],2)+ TMath::Power(Ncoll_y[k],2));
    Theta=(TMath::ATan2(Ncoll_y[k],Ncoll_x[k])); 
    Ncoll_r[k]=Ra;
    Ncoll_phi[k]=Theta;
  }
  
  if(aN_part != 0 && aN_coll !=0){
    for(int k=0;k<aN_part;k++){ 
      Nch_r[k]=Npart_r[k];
      Nch_phi[k]=Npart_phi[k];
      Nch_value[k]=((0.5)*(npp)*(1-X_hard));
    } 
    
    for(int k=0;k<aN_coll;k++){ 
      Nch_r[k+aN_part]=Ncoll_r[k];
      Nch_phi[k+aN_part]=Ncoll_phi[k];
      Nch_value[k+aN_part]=((X_hard)*npp);
    }
    
    
    for(int N=1; N<2; N++){   
      double RXA=0.0;
      double RXB=0.0;
      double RXC=0.0;
      for(int k=0; k<Total_Nch;k++){ 
	RXA=RXA+(Nch_value[k]* TMath::Power(Nch_r[k],(3.0)) );  
	RXB=RXB+(Nch_value[k]* TMath::Power(Nch_r[k],(3.0))*TMath::Cos((N)*Nch_phi[k] )  );    
	RXC=RXC+(Nch_value[k]* TMath::Power(Nch_r[k],(3.0))*TMath::Sin((N)*Nch_phi[k] )  );
      } 
      double R1=-(RXC/RXA);
      double R2=-(RXB/RXA);
      eps[N]= TMath::Sqrt((R1*R1)+(R2*R2));
      psi[N]=((TMath::ATan2(R1,R2)))/ (N);
    } //N-loop end
    	   
    for(int N = 2; N < Norder; N++){  
      double RXA=0.0;
      double RXB=0.0;
      double RXC=0.0;
      for(int k=0; k<Total_Nch;k++) { 
	RXA=RXA+(Nch_value[k]* TMath::Power(Nch_r[k],(N)));  
	RXB=RXB+(Nch_value[k]* TMath::Power(Nch_r[k],(N))*TMath::Cos((N)*Nch_phi[k] )  );    
	RXC=RXC+(Nch_value[k]* TMath::Power(Nch_r[k],(N))*TMath::Sin((N)*Nch_phi[k] )  );
      } 
      double R1=-(RXC/RXA);
      double R2=-(RXB/RXA);
      eps[N]= TMath::Sqrt((R1*R1)+(R2*R2));
      psi[N]=((TMath::ATan2(R1,R2)))/ (N);
    } 
    
  }
  else{
    for(int ii=0; ii<Norder; ii++){
      eps[ii]=0.0;
      psi[ii]=0.0;
    }
  }  //ENDIF 1
  
}



void mc_glau::get_nucleus_A(double *X1, double *Y1, double* Z1, int* PFLAG)
{
  for(int j=0;j<A;j++){
    X1[j]=XA[j]; 
    Y1[j]=YA[j];
    Z1[j]=ZA[j];
    PFLAG[j]=PROT_FLAG_A[j];
  }
}



void mc_glau::get_nucleus_B(double *X1, double *Y1, double* Z1, int* PFLAG)
{
  for(int j=0;j<B;j++){
    X1[j]=XB[j];
    Y1[j]=YB[j]; 
    Z1[j]=ZB[j]; 
    PFLAG[j]=PROT_FLAG_B[j];
  }
}





