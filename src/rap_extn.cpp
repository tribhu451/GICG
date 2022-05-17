#include "rap_extn.h"

rapidity_extension::rapidity_extension(InputData *inparams_){
iparams = inparams_ ; 
}

rapidity_extension::~rapidity_extension(){
}

void rapidity_extension::write_rapidity_extended_tilted_profile_from_idnni_boost_invariant_deposition
                 (int flag_nb, int flag_extern_file, std::string extern_file_name ){
 
  double eta_max = iparams->etamax ; 
  int neta = iparams->neta ;
  double deta = ( 2 * eta_max ) / ( neta - 1 ) ; 
  if(neta == 1){
    deta = 0 ; 
  }
  double eta_plateau = iparams->matter_eta_plateau ;
  double eta_fall = iparams->matter_eta_fall ;
  double etam = iparams->matter_tilt_param ;  

  std::string dummy_string;
  double dummy_double ; 
  int nx ;
  int ny ; 
  double dx ; 
  double dy ; 

 
  std::ifstream infile;
  if(flag_extern_file > 0){
    cout << "reading an external file for rapidity extension ..." << endl ;
    infile.open(extern_file_name.c_str());
    if (!infile){
      cout << "couldn't find external file." << endl ;
      exit(1);
    }
    else{
      cout<<"external file is : " << extern_file_name.c_str() << endl  ;
    }
  }
  else{
    cout << "reading the event averaged profile ..." << endl ;
    std::string internal_file_name = "output/idnni_boost_invariant_event_averaged_profile_for_rapidity_extension.dat" ; 
    infile.open(internal_file_name.c_str());
    if (!infile){
      cout << "couldn't find the event averaged file." << endl;
      exit(1);
    }
    else{
      cout << "the event averaged file is : " << internal_file_name.c_str() << endl  ;
    }
  }
 
   
  infile.getline(buff,200);
  iss = new istringstream(buff);
  *iss >> dummy_string >> dummy_string >> dummy_double
       >> dummy_string >> dummy_double >> dummy_string >> nx >> dummy_string >> ny
       >> dummy_string >> dummy_double >> dummy_string >> dx >> dummy_string >> dy ;
  delete iss ; 
 
  cout << "nx : " << nx << endl ; 
  cout << "ny : " << ny << endl ; 
  cout << "dx : " << dx << endl ; 
  cout << "dy : " << dy << endl ; 
 
  double x[nx][ny];
  double y[nx][ny];
  double wnet[nx][ny];
  double wa[nx][ny];
  double wb[nx][ny];

  // required while setting normalisation constant for baryon profile
  double net_participant = 0. ; 
  double participant_contribution_from_a = 0 ; 
  double participant_contribution_from_b = 0 ; 

  for(int ix = 0 ; ix < nx ; ix++ ){
    for(int iy = 0 ; iy < ny ; iy++ ){
       infile.getline(buff,200);
       iss = new istringstream(buff);
       *iss >> dummy_double >> x[ix][iy] >> y[ix][iy]
                     >> wnet[ix][iy] >> wa[ix][iy] >> wb[ix][iy] >> dummy_double >> dummy_double >> dummy_double
                     >> dummy_double  >> dummy_double  >> dummy_double  >> dummy_double ;
       net_participant += wnet[ix][iy] ;
       participant_contribution_from_a += wa[ix][iy] ; 
       participant_contribution_from_b += wb[ix][iy] ; 
       delete iss ; 
    }
  }
  infile.close();
  cout << "total participants : " << net_participant  << endl ; 

  
  cout << "extending in rapidity ... " << endl ; 
  cout << "neta : " << neta << endl ; 
  cout << "deta : " << deta << endl ; 
  if(flag_nb > 0){
    cout << "non zero baryon density inside the fluid." << endl ;  
  }

  double baryon_profile_normalisation = 0 ;
  double baryon_envelope_a = 0. ; 
  double baryon_envelope_b = 0. ;  
  double net_deposited_baryon = 0. ;

  std::ofstream outfile;
  outfile.open("output/3D_initial_profile_for_music_from_idnni.dat", std::ios::out);
  outfile<<"#"<<"\t"<<"smooth_profile_idnni"<<"\t"<<"1"<<"\t"<<"neta="<<"\t"<<neta<<"\t"<<"nx="<<"\t"<<nx<<"\t"<<"ny="<<"\t"<<ny
       <<"\t"<<"deta="<<"\t"<<deta<<"\t"<<"dx="<<"\t"<<dx<<"\t"<<"dy="<<"\t"<<dy<<endl;
  for(int ieta = 0 ; ieta < neta ; ieta++ ){
    double eta = -eta_max + deta * ieta ; 
    // cout << "ieta : " << ieta << "  eta : " << eta << endl ; 
    double plateau_profile = exp(-pow(fabs(eta) - eta_plateau, 2) / ( 2 * pow(eta_fall, 2)) * theta(fabs(eta)-eta_plateau) );
    double fa = rapidity_asymmetry_ansatz_for_tilt_by_Bozek(eta,etam);
    double fb = rapidity_asymmetry_ansatz_for_tilt_by_Bozek(-eta,etam);

    if(flag_nb > 0){
       baryon_profile_normalisation = 1.0 /  integrate_baryon_density_eta_envelop_profile_over_eta() ;
       baryon_envelope_a =  baryon_density_eta_envelop_profile(eta, 
              iparams->baryon_rapidity_profile_eta_peak, iparams->baryon_rapidity_profile_sigma_eta_plus, 
                  iparams->baryon_rapidity_profile_sigma_eta_minus ) ; 
       baryon_envelope_b = baryon_density_eta_envelop_profile(eta, 
              -iparams->baryon_rapidity_profile_eta_peak, iparams->baryon_rapidity_profile_sigma_eta_minus, 
                  iparams->baryon_rapidity_profile_sigma_eta_plus ) ; 
    }

    for(int ix = 0 ; ix < nx ; ix++ ){
      for(int iy = 0 ; iy < ny ; iy++ ){
        double matter_density =  2. * ( wa[ix][iy] * fa + wb[ix][iy] * fb ) ;
        matter_density *= plateau_profile ;

        // baryon profile
        double baryon_density  = 0. ;
        if(flag_nb > 0){
          baryon_density += baryon_envelope_a * wa[ix][iy] ;
          baryon_density += baryon_envelope_b * wb[ix][iy] ;
          baryon_density *= baryon_profile_normalisation ; 
          baryon_density /= ( dx * dy ) ;
        }
   
        double utau   = 1. ; 
        double ux     = 0. ;
        double uy     = 0. ;
        double ueta   = 0. ;
        double strange_density = 0. ;
        double charge_density  = 0. ;
        net_deposited_baryon += baryon_density ;

  
        outfile << eta << "\t" << x[ix][iy] << "\t" << y[ix][iy] << "\t" << matter_density 
                << "\t" << utau << "\t" << ux << "\t" << uy << "\t" << ueta 
                << "\t" << baryon_density << "\t" << strange_density 
                << "\t" << charge_density << endl ; 
      }
    }
  }
  outfile.close();
  cout << "net deposited baryon : " << net_deposited_baryon * dx * dy * deta << endl ; // Why have not I multiplied ( dx * dy ) ?
                                                                             // Ans : During calculation of effective-NPart, 
                                                                             // I had put all contributions to cell center
                                                                             //  in order to conserve NPart, instead of   
                                                                             // distributing it over dx*dy area around the cell.
                                                                          
                                                                             // remedy : In principle one should 
                                                                             // take the cell height to be 
                                                                             // net contribution / (dx*dy)
                                                                             // during IDNNI calculation.
  if(flag_nb > 0){ 
    if(fabs(net_participant - net_deposited_baryon * dx * dy * deta ) > 0.5 ){
       cout << "net baryon is not conserved ... " << endl ;
       exit(1); 
    }
 }
                                                                              

}



double rapidity_extension::rapidity_asymmetry_ansatz_for_tilt_by_Bozek(double etas, double etam){
  if( etas < -etam ){
   return 0.0 ;
  }
  else if( -etam <= etas && etas <= etam ){
    return ( etas + etam ) / ( 2 * etam ) ;
  }
  else{
    return 1.0;
  }
}


int rapidity_extension::theta(double xx)
{
  if(xx < 0){
    return 0;
  }
  else{
    return 1;
  }
}





///////////////////////
// Baryon deposition //
///////////////////////


double rapidity_extension::baryon_density_eta_envelop_profile(double eta, double eta0_nb, double sigma_eta_nb_plus, double sigma_eta_nb_minus ){ // as taken in arxiv:1804.10557
  double THETA_ARG = eta - eta0_nb ;
  if ( fabs(THETA_ARG) < 1E-6 ) THETA_ARG = 1E-6 ;
  return ( theta(THETA_ARG)
                    * exp( - pow( THETA_ARG , 2 ) / ( 2 * pow( sigma_eta_nb_plus, 2) ) )
                    + theta( -THETA_ARG )
                    * exp( - pow( THETA_ARG, 2 ) / ( 2 * pow( sigma_eta_nb_minus, 2) ) ) ) ;

}

double rapidity_extension::baryon_density_eta_envelop_profile_function(double* x, double* p){ // as taken in arxiv:1804.10557
  double eta_0_nb = iparams->baryon_rapidity_profile_eta_peak ;
  double sigma_eta_nb_plus = iparams->baryon_rapidity_profile_sigma_eta_plus ;
  double sigma_eta_nb_minus = iparams->baryon_rapidity_profile_sigma_eta_minus ;
  return p[0]*( theta( x[0] - eta_0_nb ) * exp( - pow( eta_0_nb - x[0], 2 ) / ( 2 * pow( sigma_eta_nb_plus, 2) ) )
                   + theta( eta_0_nb - x[0] ) * exp( - pow( eta_0_nb - x[0], 2 ) / ( 2 * pow( sigma_eta_nb_minus, 2) ) ) ) ;
}


double rapidity_extension::integrate_baryon_density_eta_envelop_profile_over_eta(){
  TF1* f1;
  f1=new TF1 ("TA",this,&rapidity_extension::baryon_density_eta_envelop_profile_function,-10,10,1);
  f1->SetParameter(0, 1.0);
  double TA = f1->Integral(-10,10,1.0E-09);
  delete f1 ;
  return TA;
}













