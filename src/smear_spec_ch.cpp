#include "smear_spec_ch.h"


smear_spec_ch::smear_spec_ch(InputData *inparams_,grid* arena_, mc_glau* _mc){
  inparams = inparams_ ; 
  mc = _mc ;
  arena = arena_ ; 
}

smear_spec_ch::~smear_spec_ch(){
}


void smear_spec_ch::smear_spectator_protons_charge(double sigma_perp, double phipp2){


  for(int ii = 0; ii < 500; ii++){
    spectator_x_of_nucleus_a[ii] = -9999. ; 
    spectator_y_of_nucleus_a[ii] =  9999. ;
    proton_flag_of_the_spectator_in_nucleus_a[ii] = 0 ;  
    spectator_x_of_nucleus_b[ii] = -9999. ; 
    spectator_y_of_nucleus_b[ii] =  9999. ;
    proton_flag_of_the_spectator_in_nucleus_b[ii] = 0 ;  
  }

 
  // collect the npart and ncoll sources info from MC Glauber class.
  int Nspec_a = mc->get_no_of_spectators_in_nucleus_a() ;  
  int Nspec_b = mc->get_no_of_spectators_in_nucleus_b() ;

  int Nspec_proton_a = 0 ; 
  int Nspec_proton_b = 0 ;  

  mc->get_spectator_positions_and_its_proton_flags_of_nucleus_a(spectator_x_of_nucleus_a,
    spectator_y_of_nucleus_a,proton_flag_of_the_spectator_in_nucleus_a ); 
  mc->get_spectator_positions_and_its_proton_flags_of_nucleus_b(spectator_x_of_nucleus_b, 
    spectator_y_of_nucleus_b,proton_flag_of_the_spectator_in_nucleus_b ); 

  for(int ii = 0; ii < Nspec_a; ii++){
    if(proton_flag_of_the_spectator_in_nucleus_a[ii] > 0){
     Nspec_proton_a += 1 ; 
    }   
  }


  for(int ii = 0; ii < Nspec_b; ii++){
    if(proton_flag_of_the_spectator_in_nucleus_b[ii] > 0){
     Nspec_proton_b += 1 ; 
    }   
  }



  // set the contribution to each cell to be zero at the beginning of the event.
  for(int ix = 0 ; ix < arena->get_nx() ; ix++ ){
    for(int iy = 0 ; iy < arena->get_ny() ; iy++ ){
      arena->get_cell(ix,iy)->set_charge_contri_from_spectator_nucleus_a_after_gaussian_smearing(0.);
      arena->get_cell(ix,iy)->set_charge_contri_from_spectator_nucleus_b_after_gaussian_smearing(0.);
    }
  }


 
  double net_ch_a = 0 ; 
  // smear the participant sources of Nucleus A.
  for(int ipart=0; ipart < Nspec_a ; ipart++ ){
    if( proton_flag_of_the_spectator_in_nucleus_a[ipart] < 1){
      continue ; 
    }
    double contributors_x ; 
    double contributors_y ; 
    contributors_x =  cos(phipp2) * spectator_x_of_nucleus_a[ipart] + sin(phipp2) * spectator_y_of_nucleus_a[ipart] ; 
    contributors_y = -sin(phipp2) * spectator_x_of_nucleus_a[ipart] + cos(phipp2) * spectator_y_of_nucleus_a[ipart] ;
    int upper_index_x, upper_index_y ; 
    get_nearest_cell_index( contributors_x + 4 * sigma_perp , contributors_y + 4 * sigma_perp , upper_index_x, upper_index_y ) ; 
    int lower_index_x, lower_index_y ;  
    get_nearest_cell_index( contributors_x - 4 * sigma_perp , contributors_y - 4 * sigma_perp , lower_index_x, lower_index_y ) ;
    for(int ix = lower_index_x ; ix < upper_index_x ; ix++ ){
      for(int iy = lower_index_y ; iy < upper_index_y ; iy++ ){
        double grid_x = -arena->get_xmax() + ix * arena->get_dx() ; 
        double grid_y = -arena->get_ymax() + iy * arena->get_dy() ;
	double weight = 0.25 * ( -erf( (grid_x - arena->get_dx() / 2 - contributors_x) / sqrt(2) / sigma_perp ) + 
				     erf((grid_x + arena->get_dx() / 2 - contributors_x) / sqrt(2) / sigma_perp) )
	    * ( -erf( (grid_y - arena->get_dy() / 2 - contributors_y) / sqrt(2) / sigma_perp ) + 
		erf((grid_y + arena->get_dy() / 2 - contributors_y) / sqrt(2) / sigma_perp) ) ;
	  
        arena->get_cell(ix,iy)->update_charge_contri_from_spectator_nucleus_a_after_gaussian_smearing(weight);
        net_ch_a += weight ; 
      } // iy loop
    } // ix loop
  } // ipart loop
  
  if(fabs(net_ch_a-Nspec_proton_a)>1){
     std::cout << "charge conservation not followed after smearing ... " << std::endl ;
     exit(-1); 
  }



   double net_ch_b = 0 ; 
  // smear the participant sources of Nucleus B.
  for(int ipart=0; ipart < Nspec_b ; ipart++ ){
    if( proton_flag_of_the_spectator_in_nucleus_b[ipart] < 1){
      continue ; 
    }
    double contributors_x ; 
    double contributors_y ; 
    contributors_x =  cos(phipp2) * spectator_x_of_nucleus_b[ipart] + sin(phipp2) * spectator_y_of_nucleus_b[ipart] ; 
    contributors_y = -sin(phipp2) * spectator_x_of_nucleus_b[ipart] + cos(phipp2) * spectator_y_of_nucleus_b[ipart] ;
    int upper_index_x, upper_index_y ; 
    get_nearest_cell_index( contributors_x + 4 * sigma_perp , contributors_y + 4 * sigma_perp , upper_index_x, upper_index_y ) ; 
    int lower_index_x, lower_index_y ;  
    get_nearest_cell_index( contributors_x - 4 * sigma_perp , contributors_y - 4 * sigma_perp , lower_index_x, lower_index_y ) ;
    for(int ix = lower_index_x ; ix < upper_index_x ; ix++ ){
      for(int iy = lower_index_y ; iy < upper_index_y ; iy++ ){
        double grid_x = -arena->get_xmax() + ix * arena->get_dx() ; 
        double grid_y = -arena->get_ymax() + iy * arena->get_dy() ;
	double weight = 0.25 * ( -erf( (grid_x - arena->get_dx() / 2 - contributors_x) / sqrt(2) / sigma_perp ) + 
				     erf((grid_x + arena->get_dx() / 2 - contributors_x) / sqrt(2) / sigma_perp) )
	    * ( -erf( (grid_y - arena->get_dy() / 2 - contributors_y) / sqrt(2) / sigma_perp ) + 
		erf((grid_y + arena->get_dy() / 2 - contributors_y) / sqrt(2) / sigma_perp) ) ;
	  
        arena->get_cell(ix,iy)->update_charge_contri_from_spectator_nucleus_b_after_gaussian_smearing(weight);
        net_ch_b += weight ; 
      } // iy loop
    } // ix loop
  } // ipart loop
  
  if(fabs(net_ch_b-Nspec_proton_b)>1){
     std::cout << "charge conservation not followed after smearing ... " << std::endl ;
     exit(-1); 
  }

  //std::cout << "Net positive charge carried by spectators of A  = " <<  net_ch_a << std::endl ; 
  //std::cout << "Net positive charge carried by spectators of B  = " <<  net_ch_b << std::endl ; 
}




void smear_spec_ch::update_charge_contribution_on_cells_over_all_events_with_spectator_gaussian_smearing(){
  for(int ix = 0 ; ix < arena->get_nx() ; ix++ ){
    for(int iy = 0 ; iy < arena->get_ny() ; iy++ ){
	arena->get_cell(ix,iy)->update_charge_contribution_from_spectator_nucleus_a_over_all_events(
            arena->get_cell(ix,iy)->get_charge_contri_from_spectator_nucleus_a_after_gaussian_smearing()) ;
	arena->get_cell(ix,iy)->update_charge_contribution_from_spectator_nucleus_b_over_all_events(
                  arena->get_cell(ix,iy)->get_charge_contri_from_spectator_nucleus_b_after_gaussian_smearing()) ;
    }
  }
}



// relevant to generate Event-by-Event profile //
void smear_spec_ch::reset_spectator_charge_contribution_from_all_events_to_zero_on_the_cells(){
  for(int ix = 0 ; ix < arena->get_nx() ; ix++ ){
    for(int iy = 0 ; iy < arena->get_ny() ; iy++ ){
      arena->get_cell(ix,iy)->reset_spectator_charge_contributions_to_zero();
    }
  }
}




// Write the event averaged profile in MUSIC format.
void smear_spec_ch::write_event_averaged_spectator_charge_profile_to_file_after_gaussian_smearing(int nEvents, int index){
  std::ofstream out_file;
  std::stringstream output_filename;
  output_filename.str("");
  output_filename << "output/event_averaged_profile_for_charge_distribution_of_spectators_" << index  ;
  output_filename << ".dat";
  out_file.open(output_filename.str().c_str(), std::ios::out);
  out_file <<"#"<<"\t"<<"event_avergaed_charge_profile"<<"\t"<<"1"<<"\t"<<"neta="<<"\t"<<"1"<<"\t"<<"nx="
            <<"\t"<<arena->get_nx()<<"\t"<<"ny="<<"\t"<<arena->get_ny()
	      <<"\t"<<"deta="<<"\t"<<"0.1"<<"\t"<<"dx="<<"\t"<<arena->get_dx()<<"\t"<<"dy="<<"\t"<<arena->get_dy()<<endl;
  double contri_of_a = 0. ; 
  double contri_of_b = 0. ; 
  double net_ch_of_spectator_from_a = 0 ;  
  double net_ch_of_spectator_from_b = 0 ;  
  for(int ix = 0 ; ix < arena->get_nx() ; ix++ ){
    for(int iy = 0 ; iy < arena->get_ny() ; iy++ ){
      double grid_x = -arena->get_xmax() + ix * arena->get_dx() ; 
      double grid_y = -arena->get_ymax() + iy * arena->get_dy() ;
        contri_of_a = arena->get_cell(ix,iy)->get_charge_contribution_from_spectator_nucleus_a_over_all_events() ; 
        contri_of_b = arena->get_cell(ix,iy)->get_charge_contribution_from_spectator_nucleus_b_over_all_events() ; 
  
        if( contri_of_a < 0.0000001 ){
	  contri_of_a = 0.0000001 ; }
        if( contri_of_b < 0.0000001 ){
	  contri_of_b = 0.0000001 ; }

      contri_of_a  /=  nEvents ; 
      contri_of_b  /=  nEvents ; 
      net_ch_of_spectator_from_a += contri_of_a ; 
      net_ch_of_spectator_from_b += contri_of_b ; 
      out_file << "0" << "  " << grid_x << "  " << grid_y << "  " << contri_of_a << "  " << contri_of_b << std::endl ;
      
    } // iy loop
  } // ix loop
  out_file.close() ; 
  std::cout << "Net-ch,  A : " <<  net_ch_of_spectator_from_a << "  B : " << net_ch_of_spectator_from_b  << std::endl ; 
}




void smear_spec_ch::write_event_averaged_spectator_charge_profile_to_file_after_gaussian_smearing_with_symmetry(int nEvents, int index){
  std::ofstream out_file;
  std::stringstream output_filename;
  output_filename.str("");
  output_filename << "output/event_averaged_profile_for_charge_distribution_of_spectators_with_symmetry_" << index  ;
  output_filename << ".dat";
  out_file.open(output_filename.str().c_str(), std::ios::out);
  
  out_file <<"#"<<"\t"<<"event_avergaed_charge_profile_w_symmetry"<<"\t"<<"1"<<"\t"<<"neta="<<"\t"<<"1"<<"\t"<<"nx="
	   <<"\t"<<arena->get_nx()<<"\t"<<"ny="<<"\t"<<arena->get_ny()
	   <<"\t"<<"deta="<<"\t"<<"0.1"<<"\t"<<"dx="<<"\t"<<arena->get_dx()<<"\t"<<"dy="<<"\t"<<arena->get_dy()<<endl;
  
  double contri_of_a = 0. ; 
  double contri_of_b = 0. ; 
  double net_ch_of_spectator_from_a = 0 ;  
  double net_ch_of_spectator_from_b = 0 ;  
  for(int ix = 0 ; ix < arena->get_nx() ; ix++ ){
    for(int iy = 0 ; iy < arena->get_ny() ; iy++ ){
      double grid_x = -arena->get_xmax() + ix * arena->get_dx() ; 
      double grid_y = -arena->get_ymax() + iy * arena->get_dy() ;
      contri_of_a = arena->get_cell(ix,iy)->get_charge_contribution_from_spectator_nucleus_a_over_all_events() ; 
      contri_of_b = arena->get_cell(arena->get_nx()-1-ix,iy)->get_charge_contribution_from_spectator_nucleus_a_over_all_events() ; 
      
      if( contri_of_a < 0.0000001 ){
	contri_of_a = 0.0000001 ; }
      if( contri_of_b < 0.0000001 ){
	contri_of_b = 0.0000001 ; }
      
      contri_of_a  /=  nEvents ; 
      contri_of_b  /=  nEvents ; 
      net_ch_of_spectator_from_a += contri_of_a ; 
      net_ch_of_spectator_from_b += contri_of_b ; 
      out_file << "0" << "  " << grid_x << "  " << grid_y << "  " << contri_of_a << "  " << contri_of_b << std::endl ;
      
    } // iy loop
  } // ix loop
  out_file.close() ; 
  std::cout << "Net-ch,  A : " <<  net_ch_of_spectator_from_a << "  B : " << net_ch_of_spectator_from_b  << std::endl ; 
}



void smear_spec_ch::write_EM_field_from_event_averaged_spectator_charge_profile_to_file_assuming_symmetry(int nEvents, int index){
  
  // Compare the strings
  if (inparams->projectile != inparams->target ) {
      std::cout << "projectile and target nuclei are not the same." 
       << "Not good to assume symmetry. Exiting the function ... "  
       << std::endl;
      return 1; 
  } else {
    std::cout << "calculating EM filed assuming symmetry" 
    <<" and storing in file ..." << std::endl ; 
  }

  double Yb = get_beam_rapidity() ;
  double tau_max = inparams->EMtaumax ;  
  double etas_max = inparams->EMetasmax ; 
  double x_max = inparams->EMxmax; 
  double y_max = inparams->EMymax; 
  
  int ntau = inparams->EMntau ; 
  int netas = inparams->EMnetas ; 
  int nxx = inparams->EMnx ; 
  int nyy = inparams->EMny ;
  
  double tau_min = 0.01 ;  
  
  double dtau = ( tau_max - tau_min ) / ( ntau ) ;
  double detas = ( 2 * etas_max ) / ( netas - 1 ) ;
  double dxx = ( 2 * x_max ) / ( nxx - 1 ) ;
  double dyy = ( 2 * y_max ) / ( nyy - 1 ) ;
  
  double emsigma = inparams->EMsigma ; 

  int MaxContributor = ( arena->get_nx() / 2 ) 
            *  ( arena->get_ny() / 2 ) ; 

  double* proj_contributor_x = new double[MaxContributor];
  double* proj_contributor_y = new double[MaxContributor];
  double* proj_contributor_charge = new double[MaxContributor];

  double* trgt_contributor_x = new double[MaxContributor];
  double* trgt_contributor_y = new double[MaxContributor];
  double* trgt_contributor_charge = new double[MaxContributor];

  for(int ii=0; ii<MaxContributor; ii++){
    proj_contributor_x[ii] = 0 ; 
    proj_contributor_y[ii] = 0 ; 
    proj_contributor_charge[ii] = 0 ; 
  }
  for(int ii=0; ii<MaxContributor; ii++){
    trgt_contributor_x[ii] = 0 ; 
    trgt_contributor_y[ii] = 0 ; 
    trgt_contributor_charge[ii] = 0 ; 
  }

   int Ncontributor_from_proj_spec = 0 ; 
   int Ncontributor_from_trgt_spec = 0 ; 
   for(int ix = 0 ; ix < arena->get_nx() ; ix++ ){
     for(int iy = 0 ; iy < arena->get_ny() ; iy++ ){
        double charge_contri = arena->get_cell(ix,iy)->
	get_charge_contribution_from_spectator_nucleus_a_over_all_events() / nEvents ; //
	if( charge_contri < 0.001 ){
          continue ; 
	}
        else{
          proj_contributor_x[Ncontributor_from_proj_spec] = -arena->get_xmax() + ix * arena->get_dx() ; 
          proj_contributor_y[Ncontributor_from_proj_spec] = -arena->get_ymax() + iy * arena->get_dy() ; 
          proj_contributor_charge[Ncontributor_from_proj_spec] = charge_contri ; 

          trgt_contributor_x[Ncontributor_from_trgt_spec] = -proj_contributor_x[Ncontributor_from_proj_spec] ; 
          trgt_contributor_y[Ncontributor_from_trgt_spec] = proj_contributor_y[Ncontributor_from_proj_spec] ; 
          trgt_contributor_charge[Ncontributor_from_trgt_spec] = charge_contri ; 

          Ncontributor_from_proj_spec += 1 ; 
          Ncontributor_from_trgt_spec += 1 ; 
        }       
     } // iy loop
   } // ix loop
  std::cout << "Ncontributor_from_projectile_spectator = " << Ncontributor_from_proj_spec << std::endl ;  
  std::cout << "Ncontributor_from_target_spectator = " << Ncontributor_from_trgt_spec << std::endl ;  

  
  std::ofstream out_file;
  std::stringstream output_filename;
  output_filename.str("");
  output_filename << "output/EM_field_profile_from_charge_distribution_" ; 
  output_filename << "of_spectators_assuming_spectators_symmetry_emsigma_" ; 
  output_filename << emsigma << "_"  << index  ;
  output_filename << ".dat";
  out_file.open(output_filename.str().c_str(), std::ios::out);
  out_file <<"#"<<"\t"<<"EM_field_profile_w_spectator_symmetry"<<"\t"<<"1"<<"\t"
	   <<"EMnetas="<<"\t"<<netas<<"\t"<<"EMnx=" <<"\t"<<nxx<<"\t"<<"EMny="
	   <<"\t"<<nyy <<"\t"<<"EMdeta="<<"\t"<< detas <<"\t"<<"EMdx="<<"\t"
	   <<dxx<<"\t"<<"EMdy="<<"\t"<<dyy<<endl;

  for(int itau = 0 ; itau < ntau ; itau++ ){
   // std::cout << "itau = " << itau << "  tau = " 
     //   << (tau_min + itau * dtau) << std::endl ; 
    for(int ietas = 0 ; ietas < netas ; ietas++ ){
      for(int ixx = 0 ; ixx < nxx ; ixx++ ){
	for(int iyy = 0 ; iyy < nyy ; iyy++ ){
	  double grid_tau = tau_min + itau * dtau ; 
	  double grid_etas = -etas_max + ietas * detas ; // eta
	  double grid_x = -x_max + ixx * dxx ; // x
	  double grid_y = -y_max + iyy * dyy ; // y
	  
	  // contribution from projectile spectators
          double By = 0 ; 
          #pragma omp parallel for reduction(+:By)
	  for(int ii = 0 ; ii < Ncontributor_from_proj_spec ; ii++ ){
	      By += get_By(proj_contributor_charge[ii],grid_tau,grid_etas,
                      grid_x,grid_y,proj_contributor_x[ii],proj_contributor_y[ii],Yb,emsigma) ; 
          }
	  // contribution from target spectators
          #pragma omp parallel for reduction(+:By)
	  for(int ii = 0 ; ii < Ncontributor_from_trgt_spec ; ii++ ){
	      By += get_By(trgt_contributor_charge[ii],grid_tau,grid_etas,
                      grid_x,grid_y,trgt_contributor_x[ii],trgt_contributor_y[ii],-Yb,emsigma) ; 
          }
	  
          if( ( fabs(grid_x)<0.001 &&  fabs(grid_y)<0.001 ) && fabs(grid_etas)<0.001 ){
            std::cout << "tau = " << grid_tau <<"  By = " << By << "  Nproj = " << Ncontributor_from_proj_spec<< std::endl ; 
	    out_file << grid_tau << "  "  << grid_etas << "  " << grid_x << "  " << grid_y << "  " << By  << std::endl ;
          }
 
	} // iyy loop
      } // ixxx loop
    } // ietas loop
  } // itau loop
  out_file.close() ; 
  

  delete proj_contributor_x ; 
  delete proj_contributor_y ; 
  delete proj_contributor_charge ; 
  delete trgt_contributor_x ; 
  delete trgt_contributor_y ; 
  delete trgt_contributor_charge ; 
  
}



double smear_spec_ch::get_beam_rapidity(){
  double yb = 0  ; 
  std::string str1 = inparams->projectile ; 
  std::string str2 = inparams->target ; 
  if(str1.compare(str2) == 0) {
    yb = acosh(inparams->SNN/(2*0.938)) ;  
  }
  else{
    std::cout << "Projectile and target are not equal."
	      << "Calculating Yb wrongly ..." << std::endl;
    exit(-1);
  }
  return yb ; 
}


double smear_spec_ch::get_By(double Q, double tau,double etas,double xx,double yy,
         double xprime,double yprime,double Yb,double emsigma){
  double xperp = sqrt(xx * xx + yy * yy);
  double phi = atan2(yy,xx);
  double xperp_prime = sqrt(xprime * xprime + yprime * yprime);
  double phi_prime = atan2(yprime,xprime);
  double Delta = get_Delta( tau, etas, xperp, phi,
               xperp_prime, phi_prime, Yb) ;
  double A = get_A(tau,etas,Yb,Delta,emsigma) ;         
  double By = ( Q / 137.) * sinh(Yb) * ( xperp * cos(phi) - xperp_prime * cos(phi_prime) ) *
              ( emsigma * fabs(sinh(Yb)) * sqrt(Delta) * 0.5 + 1 ) * exp(A) * pow(Delta,-1.5) ;        

  return By ; 
}

double smear_spec_ch::get_Delta(double tau,double etas,double xperp,double phi,
         double xperp_prime,double phi_prime,double Yb){

  double Delta = tau * tau * sinh(Yb-etas) * sinh(Yb-etas)  + xperp * xperp +
            xperp_prime * xperp_prime - 2 * xperp * xperp_prime * cos(phi-phi_prime) ; 
  return Delta ; 
}

double smear_spec_ch::get_A(double tau,double etas,double Yb,double Delta,double emsigma){
 if(Delta < 0){
  std::cout << "Delta < 0 ... " << std::endl ;
  exit(-1); 
 }
 double val = ( emsigma / 2. ) * ( tau * sinh(Yb) * sinh(Yb-etas) 
              - fabs(sinh(Yb)) * sqrt(Delta) ) ; 
 return val ; 
}










