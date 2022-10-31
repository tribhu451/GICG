#include "mc_glau_smear.h"

mc_glau_smear::mc_glau_smear(InputData *inparams_,grid* arena_, mc_glau* _mc){
  inparams = inparams_ ; 
  mc = _mc ;
  arena = arena_ ; 
}

mc_glau_smear::~mc_glau_smear(){
}

void mc_glau_smear::smear_it(double sigma_perp){

  //std::cout << "==================================================================" << std::endl ; 
  //std::cout << "No. of participants = " << mc->get_npart() << std::endl ; 
  //std::cout << "No. of participants in A = " << mc->get_no_of_participants_in_nucleus_a() << std::endl ; 
  //std::cout << "No. of participants in B = " << mc->get_no_of_participants_in_nucleus_b() << std::endl ; 


  for(int ii = 0; ii < 500; ii++){
    npart_x[ii] = -9999. ; 
    npart_y[ii] =  9999. ; 
  }

  for(int ii = 0; ii < 500; ii++){
    npart_x_of_nucleus_a[ii] = -9999. ; 
    npart_y_of_nucleus_a[ii] =  9999. ; 
  }

  for(int ii = 0; ii < 500; ii++){
    npart_x_of_nucleus_b[ii] = -9999. ; 
    npart_y_of_nucleus_b[ii] =  9999. ; 
  }

  for(int ii = 0; ii < 10000; ii++){
    ncoll_x[ii] = -9999. ; 
    ncoll_y[ii] =  9999. ; 
  }

  // collect the npart and ncoll sources info from MC Glauber class.
  int npart = mc->get_npart() ;  
  int ncoll = mc->get_ncoll() ; 
  mc->get_npart_source_positions(npart_x, npart_y); 
  mc->get_npart_source_positions_of_nucleus_a(npart_x_of_nucleus_a, npart_y_of_nucleus_a); 
  mc->get_npart_source_positions_of_nucleus_b(npart_x_of_nucleus_b, npart_y_of_nucleus_b); 
  mc->get_ncoll_source_positions(ncoll_x, ncoll_y);
  // the sources are already adjusted to give CM = (0,0)

  // print the collected info
  /*
  for(int ipart=0; ipart<npart; ipart++){
    std::cout << ipart << "  " << npart_x[ipart] << "  " << npart_y[ipart] << std::endl ; 
  }
  std::cout << "===========================================================" << std::endl ; 
  for(int icoll=0; icoll<ncoll; icoll++){
    std::cout << icoll << "  " << ncoll_x[icoll] << "  " << ncoll_y[icoll] << std::endl ; 
  }
  */
    
   // collect the value of second order participant plane (calculated ).
   int Norder = 3 ;
   double epspp[Norder];
   double phipp[Norder];
   for(int iorder = 0 ; iorder < Norder ; iorder++ ){
     epspp[iorder] = 0. ; 
     phipp[iorder] = 0. ; 
   }
   mc->calculate_eccentricity(Norder,npart,ncoll,npart_x,npart_y,ncoll_x,ncoll_y,epspp,phipp);



  // set the contribution to each cell to be zero at the beginning of the event.
  for(int ix = 0 ; ix < arena->get_nx() ; ix++ ){
    for(int iy = 0 ; iy < arena->get_ny() ; iy++ ){
      arena->get_cell(ix,iy)->set_contri_from_nucleus_a_after_gaussian_smearing(0.);
      arena->get_cell(ix,iy)->set_contri_from_nucleus_b_after_gaussian_smearing(0.);
      arena->get_cell(ix,iy)->set_contri_from_binary_collisions_after_gaussian_smearing(0.);
    }
  }


  double npart_x_for_rot_checker[500];
  double npart_y_for_rot_checker[500];
  double ncoll_x_for_rot_checker[10000];
  double ncoll_y_for_rot_checker[10000];
  for(int ii=0; ii<500; ii++){
    npart_x_for_rot_checker[ii] = 0. ; 
    npart_y_for_rot_checker[ii] = 0. ; 
  }
  for(int ii=0; ii<10000; ii++){
    ncoll_x_for_rot_checker[ii] = 0. ; 
    ncoll_y_for_rot_checker[ii] = 0. ; 
  }


  // smear the participant sources of Nucleus A.
  for(int ipart=0; ipart < mc->get_no_of_participants_in_nucleus_a() ; ipart++ ){
    double contributors_x ; 
    double contributors_y ; 
    contributors_x =  cos(phipp[2]) * npart_x_of_nucleus_a[ipart] + sin(phipp[2]) * npart_y_of_nucleus_a[ipart] ; 
    contributors_y = -sin(phipp[2]) * npart_x_of_nucleus_a[ipart] + cos(phipp[2]) * npart_y_of_nucleus_a[ipart] ;
    npart_x_for_rot_checker[ipart] = contributors_x ; 
    npart_y_for_rot_checker[ipart] = contributors_y ; 
    int upper_index_x, upper_index_y ; 
    get_nearest_cell_index( contributors_x + 4 * sigma_perp , contributors_y + 4 * sigma_perp , upper_index_x, upper_index_y ) ; 
    int lower_index_x, lower_index_y ;  
    get_nearest_cell_index( contributors_x - 4 * sigma_perp , contributors_y - 4 * sigma_perp , lower_index_x, lower_index_y ) ;
    for(int ix = lower_index_x ; ix < upper_index_x ; ix++ ){
      for(int iy = lower_index_y ; iy < upper_index_y ; iy++ ){
        double grid_x = -arena->get_xmax() + ix * arena->get_dx() ; 
        double grid_y = -arena->get_ymax() + iy * arena->get_dy() ;
	double weight_a = 0.25 * ( -erf( (grid_x - arena->get_dx() / 2 - contributors_x) / sqrt(2) / sigma_perp ) + 
				     erf((grid_x + arena->get_dx() / 2 - contributors_x) / sqrt(2) / sigma_perp) )
	    * ( -erf( (grid_y - arena->get_dy() / 2 - contributors_y) / sqrt(2) / sigma_perp ) + 
		erf((grid_y + arena->get_dy() / 2 - contributors_y) / sqrt(2) / sigma_perp) ) ;
	  
        arena->get_cell(ix,iy)->update_contri_from_nucleus_a_after_gaussian_smearing(weight_a);
      } // iy loop
    } // ix loop
  } // ipart loop


  // smear the participant sources of Nucleus B.
  for(int ipart=0; ipart < mc->get_no_of_participants_in_nucleus_b() ; ipart++ ){
    double contributors_x ; 
    double contributors_y ; 
    contributors_x =  cos(phipp[2]) * npart_x_of_nucleus_b[ipart] + sin(phipp[2]) * npart_y_of_nucleus_b[ipart] ; 
    contributors_y = -sin(phipp[2]) * npart_x_of_nucleus_b[ipart] + cos(phipp[2]) * npart_y_of_nucleus_b[ipart] ;
    npart_x_for_rot_checker[mc->get_no_of_participants_in_nucleus_a()+ipart] = contributors_x ; 
    npart_y_for_rot_checker[mc->get_no_of_participants_in_nucleus_a()+ipart] = contributors_y ; 
    int upper_index_x, upper_index_y ; 
    get_nearest_cell_index( contributors_x + 4 * sigma_perp , contributors_y + 4 * sigma_perp , upper_index_x, upper_index_y ) ; 
    int lower_index_x, lower_index_y ;  
    get_nearest_cell_index( contributors_x - 4 * sigma_perp , contributors_y - 4 * sigma_perp , lower_index_x, lower_index_y ) ;
    for(int ix = lower_index_x ; ix < upper_index_x ; ix++ ){
      for(int iy = lower_index_y ; iy < upper_index_y ; iy++ ){
        double grid_x = -arena->get_xmax() + ix * arena->get_dx() ; 
        double grid_y = -arena->get_ymax() + iy * arena->get_dy() ;
	double weight_b = 0.25 * ( -erf( (grid_x - arena->get_dx() / 2 - contributors_x) / sqrt(2) / sigma_perp ) + 
				     erf((grid_x + arena->get_dx() / 2 - contributors_x) / sqrt(2) / sigma_perp) )
	    * ( -erf( (grid_y - arena->get_dy() / 2 - contributors_y) / sqrt(2) / sigma_perp ) + 
		erf((grid_y + arena->get_dy() / 2 - contributors_y) / sqrt(2) / sigma_perp) ) ;

        arena->get_cell(ix,iy)->update_contri_from_nucleus_b_after_gaussian_smearing(weight_b);
      } // iy loop
    } // ix loop
  } // ipart loop


  // smear the binary collison sources.
  for(int icoll=0; icoll < mc->get_ncoll() ; icoll++ ){
    double contributors_x ; 
    double contributors_y ; 
    contributors_x =  cos(phipp[2]) * ncoll_x[icoll] + sin(phipp[2]) * ncoll_y[icoll] ; 
    contributors_y = -sin(phipp[2]) * ncoll_x[icoll] + cos(phipp[2]) * ncoll_y[icoll] ;
    ncoll_x_for_rot_checker[icoll] = contributors_x ; 
    ncoll_y_for_rot_checker[icoll] = contributors_y ; 
    int upper_index_x, upper_index_y ; 
    get_nearest_cell_index( contributors_x + 4 * sigma_perp , contributors_y + 4 * sigma_perp , upper_index_x, upper_index_y ) ; 
    int lower_index_x, lower_index_y ;  
    get_nearest_cell_index( contributors_x - 4 * sigma_perp , contributors_y - 4 * sigma_perp , lower_index_x, lower_index_y ) ;
    for(int ix = lower_index_x ; ix < upper_index_x ; ix++ ){
      for(int iy = lower_index_y ; iy < upper_index_y ; iy++ ){
        double grid_x = -arena->get_xmax() + ix * arena->get_dx() ; 
        double grid_y = -arena->get_ymax() + iy * arena->get_dy() ;
	double weight_ncoll = 0.25 * ( -erf( (grid_x - arena->get_dx() / 2 - contributors_x) / sqrt(2) / sigma_perp ) + 
				     erf((grid_x + arena->get_dx() / 2 - contributors_x) / sqrt(2) / sigma_perp) )
	    * ( -erf( (grid_y - arena->get_dy() / 2 - contributors_y) / sqrt(2) / sigma_perp ) + 
		erf((grid_y + arena->get_dy() / 2 - contributors_y) / sqrt(2) / sigma_perp) );
	   
        arena->get_cell(ix,iy)->update_contri_from_binary_collisions_after_gaussian_smearing(weight_ncoll);
      } // iy loop
    } // ix loop
  } // icoll loop



   // checking the rotation.
   Norder = 3 ;
   double epsppx[Norder];
   double phippx[Norder];
   for(int iorder = 0 ; iorder < Norder ; iorder++ ){
     epsppx[iorder] = 0. ; 
     phippx[iorder] = 0. ; 
   }
   mc->calculate_eccentricity(Norder,npart,ncoll,npart_x_for_rot_checker,npart_y_for_rot_checker,ncoll_x_for_rot_checker,ncoll_y_for_rot_checker,epsppx,phippx);

   if(fabs(epspp[2]-epsppx[2]) > 0.03 ){
     std::cout << "Error in Rotation(eccentricity)." << std::endl ; 
     std::cout << "before rotation, e2 = "<< epspp[2] << ",  psi2 = " << phipp[2] << std::endl ; 
     std::cout << "after rotation, e2 = " << epsppx[2] << ",   psi2 = " << phippx[2] << std::endl ; 
     exit(1);
   }
   if(fabs(phippx[2]) > 0.01 ){
     std::cout << "Error in Rotation(participant plane)." << std::endl ; 
     std::cout << "Before rotation, e2 = "<< epspp[2] << ",  psi2 = " << phipp[2] << std::endl ; 
     std::cout << "After rotation, e2 = " << epsppx[2] << ",   psi2 = " << phippx[2] << std::endl ; 
     exit(1);
   }

}



void mc_glau_smear::update_contribution_on_cells_over_all_events_with_gaussian_smearing(){
  for(int ix = 0 ; ix < arena->get_nx() ; ix++ ){
    for(int iy = 0 ; iy < arena->get_ny() ; iy++ ){
	arena->get_cell(ix,iy)->update_contribution_from_nucleus_a_over_all_events(arena->get_cell(ix,iy)->get_contri_from_nucleus_a_after_gaussian_smearing()) ;
	arena->get_cell(ix,iy)->update_contribution_from_nucleus_b_over_all_events(arena->get_cell(ix,iy)->get_contri_from_nucleus_b_after_gaussian_smearing()) ;
	arena->get_cell(ix,iy)->update_contribution_from_binary_collisions_over_all_events
                              (arena->get_cell(ix,iy)->get_contri_from_binary_collisions_after_gaussian_smearing()) ;
    }
  }
}



// Write the event averaged profile in MUSIC format.
void mc_glau_smear::write_event_averaged_profile_to_file_after_gaussian_smearing(int nEvents, int flag_to_generate_music_boost_invariant_file){
  
  std::ofstream out_file;
  if(flag_to_generate_music_boost_invariant_file > 0 ){
    out_file.open("output/mc_glauber_event_averaged_profile_for_boost_invariant_music.dat", std::ios::out);
    // cout << "the output file could be directly used in boost invariant music ... " << endl ; 
  }
  else{
    out_file.open("output/mc_glauber_boost_invariant_event_averaged_profile_for_rapidity_extension.dat", std::ios::out);
  }

  out_file <<"#"<<"\t"<<"event_avergaed_glauber"<<"\t"<<"1"<<"\t"<<"neta="<<"\t"<<"1"<<"\t"<<"nx="
            <<"\t"<<arena->get_nx()<<"\t"<<"ny="<<"\t"<<arena->get_ny()
	      <<"\t"<<"deta="<<"\t"<<"0.1"<<"\t"<<"dx="<<"\t"<<arena->get_dx()<<"\t"<<"dy="<<"\t"<<arena->get_dy()<<endl;
  
  double mult              = 0. ;
  double mult_a            = 0. ;
  double mult_b            = 0. ;
  double mult_coll         = 0. ; // energy contribution from binary collision sources.
  double participant_sum   = 0. ; 
  double participant_diff  = 0. ; 
  double npart_contri_of_a = 0. ; 
  double npart_contri_of_b = 0. ; 
  double ncoll_contri      = 0. ; 

  double total_participants_after_event_averaging = 0 ; 
  for(int ix = 0 ; ix < arena->get_nx() ; ix++ ){
    for(int iy = 0 ; iy < arena->get_ny() ; iy++ ){
      double grid_x = -arena->get_xmax() + ix * arena->get_dx() ; 
      double grid_y = -arena->get_ymax() + iy * arena->get_dy() ;
     // if( arena->get_cell(ix,iy)->get_contribution_from_nucleus_a_over_all_events() > 0.000001 && 
     //	  arena->get_cell(ix,iy)->get_contribution_from_nucleus_b_over_all_events() > 0.000001 ){
	mult = ( 0.5 * inparams->npp * ( 1 - inparams->xhard ) ) * ( arena->get_cell(ix,iy)->get_contribution_from_nucleus_a_over_all_events() +
	         arena->get_cell(ix,iy)->get_contribution_from_nucleus_b_over_all_events() ) 
                 + ( inparams->npp * ( inparams->xhard ) ) * ( arena->get_cell(ix,iy)->get_contribution_from_binary_collisions_over_all_events() ) ; 
        mult_a    = ( 0.5 * inparams->npp * (1-inparams->xhard) )   * arena->get_cell(ix,iy)->get_contribution_from_nucleus_a_over_all_events() ; 
        mult_b    = ( 0.5 * inparams->npp * (1-inparams->xhard) )   * arena->get_cell(ix,iy)->get_contribution_from_nucleus_b_over_all_events() ; 
        mult_coll = ( inparams->npp * ( inparams->xhard ) ) * arena->get_cell(ix,iy)->get_contribution_from_binary_collisions_over_all_events() ; 
        npart_contri_of_a = arena->get_cell(ix,iy)->get_contribution_from_nucleus_a_over_all_events() ; 
        npart_contri_of_b = arena->get_cell(ix,iy)->get_contribution_from_nucleus_b_over_all_events() ; 
        ncoll_contri      = arena->get_cell(ix,iy)->get_contribution_from_binary_collisions_over_all_events() ; 

      //}
      if( mult < 0.0000001 ){
	mult              = 0.00000003 ; 
	mult_a            = 0.00000001 ; 
	mult_b            = 0.00000001 ; 
	mult_coll         = 0.00000001 ; 
        npart_contri_of_a = 0.00000001 ;  
        npart_contri_of_b = 0.00000001 ;  
        ncoll_contri      = 0.00000001 ;  
      }

      mult               /=  nEvents ; 
      mult_a             /=  nEvents ; 
      mult_b             /=  nEvents ; 
      mult_coll          /=  nEvents ; 
      npart_contri_of_a  /=  nEvents ; 
      npart_contri_of_b  /=  nEvents ; 
      ncoll_contri       /=  nEvents ; 
      total_participants_after_event_averaging += (npart_contri_of_a + npart_contri_of_b) ; 

      participant_sum  += (mult_a + mult_b) ;
      participant_diff += (mult_a - mult_b) ;  
      
      if(flag_to_generate_music_boost_invariant_file  > 0 ){ 
        mult /= (arena->get_dx()*arena->get_dy()) ; 
        out_file << "0" << "\t" << grid_x << "\t" << grid_y << "\t" << mult   
                 << "\t" << "1" << "\t" << "0" << "\t" << "0" << "\t" << "0"
                 << "\t" << "0" << "\t" << "0" << "\t" << "0" << endl ;
      }
      else{
        out_file << "0" << "\t" << grid_x << "\t" << grid_y << "\t" << mult << "\t" << mult_a  << "\t" << mult_b  
                 << "\t" << mult_coll << "\t" << npart_contri_of_a << "\t" << npart_contri_of_b << "\t" <<  ncoll_contri
                 << "\t" << "1" << "\t" << "0" << "\t" << "0" << "\t" << "0"
                 << "\t" << "0" << "\t" << "0" << "\t" << "0" << endl ;
      }  
    } // iy loop
  } // ix loop

  // double ybeam = acosh(inparams->SNN / ( 2. * 0.938 ) ) ; 
  std::cout << "total participant after event averaging = " << total_participants_after_event_averaging << std::endl ; 
  std::cout << "participant assymetry : " << ( participant_diff / participant_sum )  * 100 << " %" << std::endl ; 

}

















