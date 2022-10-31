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
       baryon_profile_normalisation = 1.0 /  integrate_baryon_density_eta_envelop_profile_arxiv_1804_10557_over_eta() ;
       baryon_envelope_a =  baryon_density_eta_envelop_profile_arxiv_1804_10557(eta, 
              iparams->baryon_rapidity_profile_eta_peak, iparams->baryon_rapidity_profile_sigma_eta_plus, 
                  iparams->baryon_rapidity_profile_sigma_eta_minus ) ; 
       baryon_envelope_b = baryon_density_eta_envelop_profile_arxiv_1804_10557(eta, 
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




// rapidity Extension of the transverse profile generated from MC Glauber.
// Bozek-Wiskiel titlted fireball ansatz.
void rapidity_extension::write_rapidity_extended_tilted_profile_from_mc_glauber_boost_invariant_deposition
                 (int flag_nb, int flag_extern_file, std::string extern_file_name ){
 
  cout << "==================================" << endl ; 
  cout << "======= Rapidity Extension =======" << endl ; 
  cout << "=======       BW model     =======" << endl ; 
  cout << "==================================" << endl ; 
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
    cout << "Reading an external file for rapidity extension ..." << endl ;
    infile.open(extern_file_name.c_str());
    if (!infile){
      cout << "Couldn't find external file." << endl ;
      exit(1);
    }
    else{
      cout<<"External file is : " << extern_file_name.c_str() << endl  ;
    }
  }
  else{
    cout << "Reading the event averaged profile ..." << endl ;
    std::string internal_file_name = "output/mc_glauber_boost_invariant_event_averaged_profile_for_rapidity_extension.dat" ; 
    infile.open(internal_file_name.c_str());
    if (!infile){
      cout << "couldn't find the event averaged file." << endl;
      exit(1);
    }
    else{
      cout << "The event averaged file is : " << internal_file_name.c_str() << endl  ;
    }
  }
 
   
  infile.getline(buff,200);
  iss = new istringstream(buff);
  *iss >> dummy_string >> dummy_string >> dummy_double
       >> dummy_string >> dummy_double >> dummy_string >> nx >> dummy_string >> ny
       >> dummy_string >> dummy_double >> dummy_string >> dx >> dummy_string >> dy ;
  delete iss ; 
 
  cout << "nx   = " << nx << endl ; 
  cout << "ny   = " << ny << endl ; 
  cout << "dx   = " << dx << endl ; 
  cout << "dy   = " << dy << endl ; 

 
  double x[nx][ny];
  double y[nx][ny];
  double wnet[nx][ny];
  double wa[nx][ny];
  double wb[nx][ny];
  double wcoll[nx][ny];
  double npart_a[nx][ny]; 
  double npart_b[nx][ny]; 
  double ncoll[nx][ny];

  // required while setting normalisation constant for baryon profile
  double net_matter_contribution                      = 0 ; 
  double matter_contribution_from_a                   = 0 ; 
  double matter_contribution_from_b                   = 0 ;
  double matter_contribution_from_binary_collisions   = 0 ; 
  double npart_of_a                   = 0 ; // Npart contri from nucleus A & B, required to calculate normalisaton
  double npart_of_b                   = 0 ; // factor of baryon profile.
  double ncoll_of_ab                  = 0 ;  


  for(int ix = 0 ; ix < nx ; ix++ ){
    for(int iy = 0 ; iy < ny ; iy++ ){
       infile.getline(buff,200);
       iss = new istringstream(buff);
       *iss >> dummy_double >> x[ix][iy] >> y[ix][iy]
                     >> wnet[ix][iy] >> wa[ix][iy] >> wb[ix][iy] >> wcoll[ix][iy] 
                     >> npart_a[ix][iy] >> npart_b[ix][iy] >> ncoll[ix][iy] >> dummy_double >> dummy_double >> dummy_double
                     >> dummy_double  >> dummy_double  >> dummy_double  >> dummy_double ;
       net_matter_contribution                    += wnet[ix][iy]  ;
       matter_contribution_from_a                 += wa[ix][iy]    ; 
       matter_contribution_from_b                 += wb[ix][iy]    ; 
       matter_contribution_from_binary_collisions += wcoll[ix][iy] ; 
       npart_of_a  += npart_a[ix][iy]  ;
       npart_of_b  += npart_b[ix][iy]  ;
       ncoll_of_ab += ncoll[ix][iy] ; 
       delete iss ; 
    }
  }
  infile.close();
  // cout << "total participants : " << (npart_of_a + npart_of_b) << endl ; 

  
  cout << "Extending in rapidity ... " << endl ; 
  cout << "neta = " << neta << endl ; 
  cout << "deta = " << deta << endl ; 
  if(flag_nb > 0){
    cout << "Non zero baryon density inside the fluid." << endl ;
    if(iparams->baryon_rapidity_profile_type == 1){
      cout << "Baryon rapidty profile type is taken as in arXiv-1804.10557 ... " << endl ;
    }  
    if(iparams->baryon_rapidity_profile_type == 2){
      cout << "Baryon rapidty profile type : iiserbpr quadratic ansatz-2 ... " << endl ;
    }  
  }

  double baryon_profile_normalisation = 0. ;
  double baryon_envelope_a            = 0. ; 
  double baryon_envelope_b            = 0. ;  
  double net_deposited_baryon         = 0. ;

  if(flag_nb > 0){
    if(iparams->baryon_rapidity_profile_type == 1 ){
      double baryon_weight_temp = ( 1 - iparams->two_component_baryon_deposition_parameter_omega ) * ( npart_of_a + npart_of_b ) +
                                                      iparams->two_component_baryon_deposition_parameter_omega * 2 * ncoll_of_ab ; // 2 should be here. 
      baryon_profile_normalisation = ( npart_of_a + npart_of_b ) 
                                   / ( integrate_baryon_density_eta_envelop_profile_arxiv_1804_10557_over_eta() * baryon_weight_temp ) ;
    }
    else if(iparams->baryon_rapidity_profile_type == 2){
      double baryon_weight_temp = ( 1 - iparams->two_component_baryon_deposition_parameter_omega ) * ( npart_of_a + npart_of_b ) +
                                                      iparams->two_component_baryon_deposition_parameter_omega * 2 * ncoll_of_ab ; // 2 should be here. 
      baryon_profile_normalisation = ( npart_of_a + npart_of_b ) 
                                   / ( integrate_baryon_density_eta_envelop_profile_iiserbpr_type_2_over_eta() * baryon_weight_temp ) ;
    }
    else{
       std::cout << "Baryon rapidity profile type " << iparams->baryon_rapidity_profile_type << " doesn't exist." << std::endl ; 
       exit(1);
    }
  }

  std::ofstream outfile;
  outfile.open("output/3D_initial_profile_for_music_from_mc_glauber.dat", std::ios::out);
  outfile<<"#"<<"\t"<<"smooth_profile_idnni"<<"\t"<<"1"<<"\t"<<"neta="<<"\t"<<neta<<"\t"<<"nx="<<"\t"<<nx<<"\t"<<"ny="<<"\t"<<ny
       <<"\t"<<"deta="<<"\t"<<deta<<"\t"<<"dx="<<"\t"<<dx<<"\t"<<"dy="<<"\t"<<dy<<endl;
  for(int ieta = 0 ; ieta < neta ; ieta++ ){
    double eta = -eta_max + deta * ieta ; 
    cout << "ieta : " << ieta << "  eta : " << eta << endl ; 
    double plateau_profile = exp(-pow(fabs(eta) - eta_plateau, 2) / ( 2 * pow(eta_fall, 2)) * theta(fabs(eta)-eta_plateau) );
    double fa = rapidity_asymmetry_ansatz_for_tilt_by_Bozek(eta,etam);
    double fb = rapidity_asymmetry_ansatz_for_tilt_by_Bozek(-eta,etam);

    if(flag_nb > 0 && iparams->baryon_rapidity_profile_type == 1 ){
       baryon_envelope_a =  baryon_density_eta_envelop_profile_arxiv_1804_10557(eta, 
              iparams->baryon_rapidity_profile_eta_peak, iparams->baryon_rapidity_profile_sigma_eta_plus, 
                  iparams->baryon_rapidity_profile_sigma_eta_minus ) ; 
       baryon_envelope_b = baryon_density_eta_envelop_profile_arxiv_1804_10557(eta, 
              -iparams->baryon_rapidity_profile_eta_peak, iparams->baryon_rapidity_profile_sigma_eta_minus, 
                  iparams->baryon_rapidity_profile_sigma_eta_plus ) ; 
    }

  
    if(flag_nb > 0 && iparams->baryon_rapidity_profile_type == 2 ){
       baryon_envelope_a =  baryon_density_eta_envelop_profile_iiserbpr_type_2(eta, 
              iparams->iiserbpr_ra_type_2_peak, iparams->iiserbpr_ra_type_2_etal, 
                  iparams->iiserbpr_ra_type_2_right_fall ) ; 
       baryon_envelope_b = baryon_density_eta_envelop_profile_iiserbpr_type_2(-eta, 
               iparams->iiserbpr_ra_type_2_peak, iparams->iiserbpr_ra_type_2_etal, 
                  iparams->iiserbpr_ra_type_2_right_fall ) ; 
    }
 


    for(int ix = 0 ; ix < nx ; ix++ ){
      for(int iy = 0 ; iy < ny ; iy++ ){
        double matter_density =  2. * ( wa[ix][iy] * fa + wb[ix][iy] * fb ) + wcoll[ix][iy]  ;
        matter_density *= plateau_profile ;
        matter_density /= ( dx * dy ) ;  // dividing by a constant (dx*dy). During smearing (dx*dy) is multiplied which is scaled down now.

        // baryon profile
        double baryon_density  = 0. ;
        if(flag_nb > 0 ){
          baryon_density += baryon_envelope_a * ( 1 - iparams->two_component_baryon_deposition_parameter_omega ) * npart_a[ix][iy] ;
          baryon_density += baryon_envelope_b * ( 1 - iparams->two_component_baryon_deposition_parameter_omega ) * npart_b[ix][iy] ;
          baryon_density += ( baryon_envelope_a + baryon_envelope_b ) * iparams->two_component_baryon_deposition_parameter_omega * ncoll[ix][iy] ;
          baryon_density *= baryon_profile_normalisation ; 
          baryon_density /= ( dx * dy ) ;
        }
   
        double utau   = 1. ; 
        double ux     = 0. ;
        double uy     = 0. ;
        double ueta   = 0. ;
        double strange_density = 0. ;
        double charge_density  = 0. ;
        net_deposited_baryon += baryon_density * deta * dx * dy  ;

  
        outfile << eta << "\t" << x[ix][iy] << "\t" << y[ix][iy] << "\t" << matter_density 
                << "\t" << utau << "\t" << ux << "\t" << uy << "\t" << ueta 
                << "\t" << baryon_density << "\t" << strange_density 
                << "\t" << charge_density << endl ; 
      }
    }
  }
  outfile.close();
  cout << "Net deposited baryon = " << net_deposited_baryon  << endl ; 

  if(flag_nb > 0){ 
    if(fabs(npart_of_a + npart_of_b - net_deposited_baryon  ) > 2.0 ){
       std::cout << "Net baryon is not conserved ... " << std::endl ;
       std::cout << "npart of a = " << npart_of_a << ",  npart of b = " << npart_of_b << ",  sum of baryon density = " << net_deposited_baryon << std::endl ; 
       exit(1); 
    }
 }
                                                                              
}







// rapidity Extension of the transverse profile generated from MC Glauber.
// SARJ = Shen-Alzhrani-Ryu-Jupic = 'Chun Shen', 'Shar Alzhrani', 'Sangwook Ryu', 'Vahidin Jupic'
// Refer to : arXiv:2106.08125 and arXiv:2003.05852
void rapidity_extension::write_rapidity_extended_energy_momentum_conserving_SARJ_profile_from_mc_glauber_boost_invariant_deposition
                 (int flag_nb, int flag_extern_file, std::string extern_file_name ){
 
  cout << "==================================" << endl ; 
  cout << "======= Rapidity Extension =======" << endl ; 
  cout << "=======     SARJ model     =======" << endl ; 
  cout << "==================================" << endl ; 
  double eta_max = iparams->etamax ; 
  int neta = iparams->neta ;
  double deta = ( 2 * eta_max ) / ( neta - 1 ) ; 
  if(neta == 1){
    deta = 0 ; 
  }
  double eta_plateau = iparams->matter_eta_plateau ;
  double eta_fall    = iparams->matter_eta_fall ;
  double f           = iparams->SARJ_model_f_param ; 
  double mN          = 0.938 ;                                    // mass of nucleon in GeV
  double ybeam       = acosh( iparams->SNN / ( 2. * mN ) ) ;      // beam rapidity. s_{NN} in GeV .  
  double ceta = 0 ; 

  std::string dummy_string;
  double dummy_double ; 
  int    nx ;
  int    ny ; 
  double dx ; 
  double dy ; 

 
  std::ifstream infile;
  if(flag_extern_file > 0){
    cout << "Reading an external file for rapidity extension ..." << endl ;
    infile.open(extern_file_name.c_str());
    if (!infile){
      cout << "Couldn't find external file." << endl ;
      exit(1);
    }
    else{
      cout<<"External file is : " << extern_file_name.c_str() << endl  ;
    }
  }
  else{
    cout << "Reading the event averaged profile ..." << endl ;
    std::string internal_file_name = "output/mc_glauber_boost_invariant_event_averaged_profile_for_rapidity_extension.dat" ; 
    infile.open(internal_file_name.c_str());
    if (!infile){
      cout << "couldn't find the event averaged file." << endl;
      exit(1);
    }
    else{
      cout << "The event averaged file is : " << internal_file_name.c_str() << endl  ;
    }
  }
 
   
  infile.getline(buff,200);
  iss = new istringstream(buff);
  *iss >> dummy_string >> dummy_string >> dummy_double
       >> dummy_string >> dummy_double >> dummy_string >> nx >> dummy_string >> ny
       >> dummy_string >> dummy_double >> dummy_string >> dx >> dummy_string >> dy ;
  delete iss ; 
 
  cout << "nx   = " << nx << endl ; 
  cout << "ny   = " << ny << endl ; 
  cout << "dx   = " << dx << endl ; 
  cout << "dy   = " << dy << endl ; 

 
  double x         [nx][ny];
  double y         [nx][ny];
  double npart_a   [nx][ny]; 
  double npart_b   [nx][ny]; 
  double M         [nx][ny];
  double ycm       [nx][ny];

  // required while setting normalisation constant for baryon profile
  double npart_of_a                                   = 0 ; // Npart contri from nucleus A & B, required to calculate normalisaton
  double npart_of_b                                   = 0 ; // factor of baryon profile.

  for(int ix = 0 ; ix < nx ; ix++ ){
    for(int iy = 0 ; iy < ny ; iy++ ){
       infile.getline(buff,200);
       iss = new istringstream(buff);
       *iss >> dummy_double >> x[ix][iy] >> y[ix][iy]
                     >> dummy_double >> dummy_double >> dummy_double >> dummy_double 
                     >> npart_a[ix][iy] >> npart_b[ix][iy] >> dummy_double >> dummy_double >> dummy_double >> dummy_double
                     >> dummy_double  >> dummy_double  >> dummy_double  >> dummy_double ;
       npart_of_a                                   +=   npart_a [ix][iy] ;
       npart_of_b                                   +=   npart_b [ix][iy] ;
       M[ix][iy]    =  mN * sqrt( npart_a[ix][iy] * npart_a[ix][iy] +  npart_b[ix][iy] * npart_b[ix][iy] +
                                   2 * npart_a[ix][iy] * npart_b[ix][iy] * cosh(2*ybeam) ) ;                                 // M(x,y) = invariant mass.
       ycm[ix][iy]  =  atanh( ( npart_a[ix][iy] - npart_b[ix][iy] ) / ( npart_a[ix][iy] + npart_b[ix][iy] ) * tanh(ybeam) ); // ycm    = center of mass rapidity.
       ceta         =  exp( eta_plateau ) * TMath::Erfc( -eta_fall / sqrt(2.) ) + 
                                   exp( -eta_plateau ) * TMath::Erfc( eta_fall / sqrt(2.) ) ;

       M[ix][iy]    /=  iparams->init_time_tau0 *
                       ( 2 * sinh(eta_plateau) + sqrt( TMath::Pi() / 2. ) * eta_fall * exp( eta_fall * eta_fall / 2 ) * ceta ) ;

       delete iss ; 
    }
  }
  infile.close();
  // cout << "total participants : " << (npart_of_a + npart_of_b) << endl ; 

  
  cout << "Extending in rapidity ... " << endl ; 
  cout << "neta = " << neta << endl ; 
  cout << "deta = " << deta << endl ; 
  if(flag_nb > 0){
    cout << "Non zero baryon density inside the fluid." << endl ;
    if(iparams->baryon_rapidity_profile_type == 1){
      cout << "Baryon rapidty profile type is taken as in arXiv-1804.10557 ... " << endl ;
    }  
    if(iparams->baryon_rapidity_profile_type == 2){
      cout << "Baryon rapidty profile type : iiserbpr quadratic ansatz-2 ... " << endl ;
    }  
  }

  double baryon_profile_normalisation = 0. ;
  double baryon_envelope_a            = 0. ; 
  double baryon_envelope_b            = 0. ;  
  double net_deposited_baryon         = 0. ;

  if(iparams->baryon_rapidity_profile_type == 1 ){
    baryon_profile_normalisation = 1. / ( integrate_baryon_density_eta_envelop_profile_arxiv_1804_10557_over_eta()  ) ;
  }
  else if(iparams->baryon_rapidity_profile_type == 2){
    baryon_profile_normalisation = 1. / ( integrate_baryon_density_eta_envelop_profile_arxiv_1804_10557_over_eta()  ) ;
  }
  else{
     std::cout << "Baryon rapidity profile type " << iparams->baryon_rapidity_profile_type << " doesn't exist." << std::endl ; 
     exit(1);
  }


  std::ofstream outfile;
  outfile.open("output/3D_initial_profile_for_music_from_mc_glauber.dat", std::ios::out);
  outfile<<"#"<<"\t"<<"smooth_profile_idnni"<<"\t"<<"1"<<"\t"<<"neta="<<"\t"<<neta<<"\t"<<"nx="<<"\t"<<nx<<"\t"<<"ny="<<"\t"<<ny
       <<"\t"<<"deta="<<"\t"<<deta<<"\t"<<"dx="<<"\t"<<dx<<"\t"<<"dy="<<"\t"<<dy<<endl;
  for(int ieta = 0 ; ieta < neta ; ieta++ ){
    double eta = -eta_max + deta * ieta ; 
    cout << "ieta : " << ieta << "  eta : " << eta << endl ; 



    if(flag_nb > 0 && iparams->baryon_rapidity_profile_type == 1 ){
       baryon_envelope_a =  baryon_density_eta_envelop_profile_arxiv_1804_10557(eta, 
              iparams->baryon_rapidity_profile_eta_peak, iparams->baryon_rapidity_profile_sigma_eta_plus, 
                  iparams->baryon_rapidity_profile_sigma_eta_minus ) ; 
       baryon_envelope_b = baryon_density_eta_envelop_profile_arxiv_1804_10557(eta, 
              -iparams->baryon_rapidity_profile_eta_peak, iparams->baryon_rapidity_profile_sigma_eta_minus, 
                  iparams->baryon_rapidity_profile_sigma_eta_plus ) ; 
    }

  
    if(flag_nb > 0 && iparams->baryon_rapidity_profile_type == 2 ){
       baryon_envelope_a =  baryon_density_eta_envelop_profile_iiserbpr_type_2(eta, 
              iparams->iiserbpr_ra_type_2_peak, iparams->iiserbpr_ra_type_2_etal, 
                  iparams->iiserbpr_ra_type_2_right_fall ) ; 
       baryon_envelope_b = baryon_density_eta_envelop_profile_iiserbpr_type_2(-eta, 
               iparams->iiserbpr_ra_type_2_peak, iparams->iiserbpr_ra_type_2_etal, 
                  iparams->iiserbpr_ra_type_2_right_fall ) ; 
    }
 


    for(int ix = 0 ; ix < nx ; ix++ ){
      for(int iy = 0 ; iy < ny ; iy++ ){
        double plateau_profile = exp(-pow(fabs(eta - (1-f)*ycm[ix][iy] ) - eta_plateau, 2) / ( 2 * pow(eta_fall, 2)) *
                                     theta(fabs(eta - (1-f)*ycm[ix][iy] ) - eta_plateau) );
        double matter_density = M[ix][iy]  ;
        matter_density *= plateau_profile ;
        matter_density /= ( dx * dy ) ;  // dividing by a constant (dx*dy). During smearing (dx*dy) is multiplied which is scaled down now.

        // baryon profile
        double baryon_density  = 0. ;
        if(flag_nb > 0 ){
          baryon_density += baryon_envelope_a * ( 1 - iparams->two_component_baryon_deposition_parameter_omega ) * npart_a[ix][iy] ;
          baryon_density += baryon_envelope_b * ( 1 - iparams->two_component_baryon_deposition_parameter_omega ) * npart_b[ix][iy] ;
          baryon_density *= baryon_profile_normalisation ; 
          baryon_density /= ( dx * dy ) ;
        }
   
        double strange_density = 0. ;
        double charge_density  = 0. ;
        net_deposited_baryon += baryon_density * deta * dx * dy  ;

        // utau, ux,uy will be set inside MUSIC code.
        outfile << eta << "\t" << x[ix][iy] << "\t" << y[ix][iy] << "\t" << matter_density 
                << "\t" << f*ycm[ix][iy] << "\t" << "0" << "\t" << "0" << "\t" << "0" 
                << "\t" << baryon_density << "\t" << strange_density 
                << "\t" << charge_density << endl ; 
      }
    }
  }
  outfile.close();
  cout << "Net deposited baryon = " << net_deposited_baryon  << endl ; 

  if(flag_nb > 0){ 
    if(fabs(npart_of_a + npart_of_b - net_deposited_baryon  ) > 2.0 ){
       std::cout << "Net baryon is not conserved ... " << std::endl ;
       std::cout << "npart of a = " << npart_of_a << ",  npart of b = " << npart_of_b << ",  sum of baryon density = " << net_deposited_baryon << std::endl ; 
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
  if(xx < 0 ){
    return 0;
  }
  else{
    return 1;
  }
}





///////////////////////
// Baryon deposition //
///////////////////////


double rapidity_extension::baryon_density_eta_envelop_profile_arxiv_1804_10557(double eta, double eta0_nb, double sigma_eta_nb_plus, double sigma_eta_nb_minus ){ 
// as taken in arxiv:1804.10557
  double THETA_ARG = eta - eta0_nb ;
  if ( fabs(THETA_ARG) < 1E-6 ) THETA_ARG = 1E-6 ;
  return ( theta(THETA_ARG)
                    * exp( - pow( THETA_ARG , 2 ) / ( 2 * pow( sigma_eta_nb_plus, 2) ) )
                    + theta( -THETA_ARG )
                    * exp( - pow( THETA_ARG, 2 ) / ( 2 * pow( sigma_eta_nb_minus, 2) ) ) ) ;

}

double rapidity_extension::baryon_density_eta_envelop_profile_arxiv_1804_10557_function(double* x, double* p){ 
// as taken in arxiv:1804.10557
  double eta_0_nb = iparams->baryon_rapidity_profile_eta_peak ;
  double sigma_eta_nb_plus = iparams->baryon_rapidity_profile_sigma_eta_plus ;
  double sigma_eta_nb_minus = iparams->baryon_rapidity_profile_sigma_eta_minus ;
  return p[0]*( theta( x[0] - eta_0_nb ) * exp( - pow( eta_0_nb - x[0], 2 ) / ( 2 * pow( sigma_eta_nb_plus, 2) ) )
                   + theta( eta_0_nb - x[0] ) * exp( - pow( eta_0_nb - x[0], 2 ) / ( 2 * pow( sigma_eta_nb_minus, 2) ) ) ) ;
}


double rapidity_extension::integrate_baryon_density_eta_envelop_profile_arxiv_1804_10557_over_eta(){
  TF1* f1;
  f1=new TF1 ("TA",this,&rapidity_extension::baryon_density_eta_envelop_profile_arxiv_1804_10557_function,-10,10,1);
  f1->SetParameter(0, 1.0);
  double TA = f1->Integral(-10,10,1.0E-06);
  delete f1 ;
  return TA;
}


// ============================================================================
// =============== IISER Berhampur baryon rapidity ansatz-2 ===================
// ============================================================================
double rapidity_extension::baryon_density_eta_envelop_profile_iiserbpr_type_2(double eta, double eta_0, double eta_L, double sig){
     if (eta > 0 ) eta_0 -= 0.0001 ;
     if (eta < 0 ) eta_0 += 0.0001 ;  
     return theta( eta - eta_0 ) * exp( -pow(eta-eta_0,2) / (2*sig*sig) ) +
       theta(eta_0 - eta) * theta(eta-eta_L) * pow( (eta-eta_L)/(eta_0-eta_L), 2 )  ;
}


double rapidity_extension::baryon_density_eta_envelop_profile_iiserbpr_type_2_function(double* x, double* p){
   double eta_0 = iparams->iiserbpr_ra_type_2_peak ; 
   if(x[0] > 0 ) eta_0 = iparams->iiserbpr_ra_type_2_peak - 0.0001 ;
   if(x[0] < 0 ) eta_0 = iparams->iiserbpr_ra_type_2_peak + 0.0001 ;
   double eta_L = iparams->iiserbpr_ra_type_2_etal ;
   double sig   = iparams->iiserbpr_ra_type_2_right_fall  ;
   return p[0] * theta( x[0] - eta_0 ) * exp( -pow(x[0]-eta_0,2) / (2*sig*sig) ) +
       theta(eta_0 - x[0]) * theta(x[0]-eta_L) * pow ( (x[0]-eta_L)/(eta_0-eta_L), 2 )   ;

}


double rapidity_extension::integrate_baryon_density_eta_envelop_profile_iiserbpr_type_2_over_eta(){
  TF1* f2;
  f2=new TF1 ("TB",this,&rapidity_extension::baryon_density_eta_envelop_profile_iiserbpr_type_2_function,-10,10,1);
  f2->SetParameter(0, 1.0);
  double TB = f2->Integral(-10,10,1.0E-06);
  delete f2;
  return TB;
}













