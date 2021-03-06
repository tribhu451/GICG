#include "job.h"


job::job(InputData* iparams_, grid* arena_, mc_glau* MC_ ){
  iparams = iparams_ ; 
  arena = arena_ ; 
  MC = MC_ ;
  if(iparams->mode == 0 ){ 
    IDNNI = new idnni(iparams,arena,MC);
  }
  if(iparams->mode == 1){
    MCGlbS = new mc_glau_smear(iparams,arena,MC);
  }
}


job::~job(){
}


void job::event_by_event_mc_glauber_for_centrality_determination(int event_no){

  string collision_species = iparams->projectile +"_"+ iparams->target +"_" + to_string(int(iparams->SNN)) ;
  double npart,ncoll,nch,b ;
  TFile* root_file = new TFile (iparams->root_output_file_name.c_str(),"recreate");
  TTree* root_tree= new TTree(collision_species.c_str(),collision_species.c_str());
  root_tree->Branch("npart",&npart,"npart/D");
  root_tree->Branch("ncoll",&ncoll,"ncoll/D");
  root_tree->Branch("nch",&nch,"nch/D");
  root_tree->Branch("b",&b,"b/D");
  
  int event_count = 0 ;
  // event loop starts
  do 
    {
      // generate the MC Glauber event. 
      MC->event(0);

      // get the values.
      npart =  MC->get_npart();
      ncoll =  MC->get_ncoll();
      nch   =  iparams->npp * ( (1.0-iparams->xhard)*0.5*MC->get_npart() 
      		      + iparams->xhard*MC->get_ncoll() ) ;
      b = MC->get_impactf();

      // if atleast two nparts are there then count that as an event.
      if (npart >= 2)
	{
	  cout << "event count : "  << event_count << endl ; 
          cout << "b : " << b << ", NPart : " << npart << ", NCh : " <<  nch << " ..."<< endl ; 
	  root_tree->Fill();
	  event_count++;
	}
    } 
  while( event_count < event_no ) ;
    
  root_tree->Write();
  root_file->Close(); 
}



void job::event_by_event_mc_glauber_and_eccentricity_calculation(int event_no){
  string collision_species = iparams->projectile +"_"+ iparams->target +"_" + to_string(int(iparams->SNN)) ;
  // root file init
  double npart,ncoll,nch,b, e1, phi1, e2, phi2, e3, phi3, e4, phi4, e5, phi5, e6, phi6 ;
  TFile* root_file = new TFile (iparams->root_output_file_name.c_str(),"recreate");
  TTree* root_tree= new TTree(collision_species.c_str(),collision_species.c_str());
  root_tree->Branch("npart",&npart,"npart/D");
  root_tree->Branch("ncoll",&ncoll,"ncoll/D");
  root_tree->Branch("nch",&nch,"nch/D");
  root_tree->Branch("b",&b,"b/D");
  root_tree->Branch("e1",&e1,"e1/D");
  root_tree->Branch("phi1",&phi1,"phi1/D");
  root_tree->Branch("e2",&e2,"e2/D");
  root_tree->Branch("phi2",&phi2,"phi2/D");
  root_tree->Branch("e3",&e3,"e3/D");
  root_tree->Branch("phi3",&phi3,"phi3/D");
  root_tree->Branch("e4",&e4,"e4/D");
  root_tree->Branch("phi4",&phi4,"phi4/D");
  root_tree->Branch("e5",&e5,"e5/D");
  root_tree->Branch("phi5",&phi5,"phi5/D");
  root_tree->Branch("e6",&e6,"e6/D");
  root_tree->Branch("phi6",&phi6,"phi6/D");
  
  int event_count = 0 ;
  // event loop starts
  do 
    {
      // generate the MC Glauber event. 
      MC->event(1);

      // get the values.
      npart =  MC->get_npart();
      ncoll =  MC->get_ncoll();
      nch   =  iparams->npp * ( (1.0-iparams->xhard)*0.5*MC->get_npart() 
      		      + iparams->xhard*MC->get_ncoll() ) ;
      b = MC->get_impactf();
      e1 = MC->eccen1();
      phi1 = MC->phi1();
      e2 = MC->eccen2();
      phi2 = MC->phi2();
      e3 = MC->eccen3();
      phi3 = MC->phi3();
      e4 = MC->eccen4();
      phi4 = MC->phi4();
      e5 = MC->eccen5();
      phi5 = MC->phi5();
      e6 = MC->eccen6();
      phi6 = MC->phi6();

      // if atleast two nparts are there then count that as an event.
      if (npart >= 2)
	{
	  cout << "event count : "  << event_count << endl ; 
          cout << "b : " << b << ", NPart : " << npart << ", NCh : " <<  nch << " ..."<< endl ; 
	  root_tree->Fill();
	  event_count++;
	}
    } 
  while( event_count < event_no ) ;
    
  root_tree->Write();
  root_file->Close(); 
  
}



void job::event_by_event_idnni_for_centrality_determination(int event_no){
  string collision_species = iparams->projectile +"_"+ iparams->target +"_" + to_string(int(iparams->SNN)) ;
  // root file init
  double nch,b ;
  TFile* root_file = new TFile (iparams->root_output_file_name.c_str(),"recreate");
  TTree* root_tree= new TTree(collision_species.c_str(),collision_species.c_str());
  root_tree->Branch("nch",&nch,"nch/D");
  root_tree->Branch("b",&b,"b/D");

  int event_count = 0 ;
  do 
    {
      // generate the MC Glauber event. 
      // It is must in all cases to generate nucleons.
      MC->event(0);
      b = MC->get_impactf();
      IDNNI->perform();

      if( IDNNI->get_total_deposition_over_all_cells_for_current_event() > 0.1 )
	{
	  cout << "event count : "  << event_count << endl ; 
          nch = IDNNI->get_total_deposition_over_all_cells_for_current_event() ;
          cout << "b : " << b << ", effective npart : " <<  nch  << " ..."<< endl ; 
	  root_tree->Fill();
	  event_count++;
	}
    } 
  while( event_count < event_no ) ;
  root_tree->Write();
  root_file->Close(); 
  
}


void job::event_by_event_idnni_and_eccentricity_calculation(int event_no){

  string collision_species = iparams->projectile +"_"+ iparams->target +"_" + to_string(int(iparams->SNN)) ;
  // root file init
  double nch,b, e1, phi1, e2, phi2, e3, phi3, e4, phi4, e5, phi5, e6, phi6 ;
  TFile* root_file = new TFile (iparams->root_output_file_name.c_str(),"recreate");
  TTree* root_tree= new TTree(collision_species.c_str(),collision_species.c_str());
  root_tree->Branch("nch",&nch,"nch/D");
  root_tree->Branch("b",&b,"b/D");
  root_tree->Branch("e1",&e1,"e1/D");
  root_tree->Branch("phi1",&phi1,"phi1/D");
  root_tree->Branch("e2",&e2,"e2/D");
  root_tree->Branch("phi2",&phi2,"phi2/D");
  root_tree->Branch("e3",&e3,"e3/D");
  root_tree->Branch("phi3",&phi3,"phi3/D");
  root_tree->Branch("e4",&e4,"e4/D");
  root_tree->Branch("phi4",&phi4,"phi4/D");
  root_tree->Branch("e5",&e5,"e5/D");
  root_tree->Branch("phi5",&phi5,"phi5/D");
  root_tree->Branch("e6",&e6,"e6/D");
  root_tree->Branch("phi6",&phi6,"phi6/D");

  int event_count = 0 ;

  int Norder = 7 ;
  double eps[Norder];
  double phi[Norder];
  for(int iorder = 0 ; iorder < Norder ; iorder++ ){
    eps[iorder] = 0. ; 
    phi[iorder] = 0. ; 
  }

  do 
    {
      // generate the MC Glauber event. 
      // It is must in all cases to generate nucleons.
      MC->event(0);
      b = MC->get_impactf();
      IDNNI->perform();

      if( IDNNI->get_total_deposition_over_all_cells_for_current_event() > 0.1 )
	{
	  cout << "event count : "  << event_count << endl ; 
          nch = IDNNI->get_total_deposition_over_all_cells_for_current_event() ;
          IDNNI->get_eccentricities_and_participant_plane_angles(Norder, eps, phi) ;
          e1 = eps[1] ; 
          e2 = eps[2] ;
          e3 = eps[3] ;
          e4 = eps[4] ;
          e5 = eps[5] ;
          e6 = eps[6] ; 
          phi1 = phi[1] ; 
          phi2 = phi[2] ; 
          phi3 = phi[3] ; 
          phi4 = phi[4] ; 
          phi5 = phi[5] ; 
          phi6 = phi[6] ; 
          cout << "b : " << b << ", effective npart : " <<  nch << ",  e2 : " << e2 << ",  psi2 : " << phi2 << " ..."<< endl ; 
	  root_tree->Fill();
	  event_count++;
	}
    } 
  while( event_count < event_no ) ;
    
  root_tree->Write();
  root_file->Close(); 
}


void job::rotate_by_second_order_participant_plane_and_then_gaussian_smearing_for_event_averaged_profile_using_idnni(int event_no){

  double b;
  double e2 ; 
  double phi2 ;
  double eps2 ; 
  double psi2 ; 

  int Norder = 3 ;
  double eps[Norder];
  double phi[Norder];
  for(int iorder = 0 ; iorder < Norder ; iorder++ ){
    eps[iorder] = 0. ; 
    phi[iorder] = 0. ; 
  }

  int event_count = 0 ;
  do 
    {
      // It is must in all cases to generate nucleons.
      MC->event(0);
      b = MC->get_impactf();
      IDNNI->perform();
      if( IDNNI->get_total_deposition_over_all_cells_for_current_event() > iparams->lower_mult_proxy_cut && 
           IDNNI->get_total_deposition_over_all_cells_for_current_event() <= iparams->upper_mult_proxy_cut )
	{
	  cout << "event count : "  << event_count << endl ; 
          IDNNI->get_eccentricities_and_participant_plane_angles(Norder, eps, phi) ;
          e2 = eps[2] ;
          phi2 = phi[2] ; 
          cout << "b : " << b << ", Net deposition : " << IDNNI->get_total_deposition_over_all_cells_for_current_event()
                       << ",  e2 : " << e2 << ",  psi2 : " << phi2 << endl ;

          IDNNI->rotation_by_second_order_event_plane_followed_by_gaussian_smearing_for_the_current_event(iparams->gaussian_smearing_sigma, phi2);
          IDNNI->get_eccentricities_and_participant_plane_angles_after_gaussian_smearing(Norder, eps, phi) ;
          eps2 = eps[2] ;
          psi2 = phi[2] ; 
          cout << "  after smear : " <<  IDNNI->get_total_deposition_over_all_cells_after_gaussian_smearing_for_current_event() 
                             << ",  e2 : " << eps2 << ",  psi2 : " << psi2 << " ..."<< endl ;
          
          if( fabs(e2-eps2) > 0.01 || fabs(psi2)> 0.001 ){
           cout << "[Warning] Wrong result after gaussian smearing ..." << endl ; 
           // exit(1);
          }
          
	  event_count++;
          IDNNI->update_contribution_on_cells_over_all_events_with_gaussian_smearing();
	}
    } 
  while( event_count < event_no ) ;
  IDNNI->write_event_averaged_profile_to_file_after_gaussian_smearing(event_no,0);
  IDNNI->write_event_averaged_profile_to_file_after_gaussian_smearing(event_no,1);

}



void job::rotate_by_second_order_participant_plane_and_then_gaussian_smearing_for_event_averaged_profile_using_mc_glauber(int event_no){

  double b;
  int event_count = 0 ;
  do 
    {
      // It is must in all cases to generate nucleons.
      MC->event(0);
      b = MC->get_impactf();
      if( MC->get_two_component_galuber_multiplicity_proxy() > iparams->lower_mult_proxy_cut && 
             MC->get_two_component_galuber_multiplicity_proxy() <= iparams->upper_mult_proxy_cut ){

           if( (event_count%50) == 0 ){
             std::cout << "event no = " << event_count << ",  b = " << b << " (fm),  Npart = " << MC->get_npart() << std::endl ; 
           }
           MCGlbS->smear_it(iparams->gaussian_smearing_sigma);
           MCGlbS->update_contribution_on_cells_over_all_events_with_gaussian_smearing();
           event_count++;
      }
    }
  while( event_count < event_no ) ;

  // [About below 2 lines] The second argument is a flag for making a boost invariant profile for music input
  MCGlbS->write_event_averaged_profile_to_file_after_gaussian_smearing(event_no,0); 
  MCGlbS->write_event_averaged_profile_to_file_after_gaussian_smearing(event_no,1);

}









