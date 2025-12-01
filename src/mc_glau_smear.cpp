#include "mc_glau_smear.h"

mc_glau_smear::mc_glau_smear(InputData *inparams_,grid* arena_, mc_glau* _mc){
  inparams = inparams_ ; 
  mc = _mc ;
  arena = arena_ ; 
  TINY = 1e-10 ; 
}

mc_glau_smear::~mc_glau_smear(){
}

void mc_glau_smear::smear_it(double sigma_perp){

  for(int ii = 0; ii < 500; ii++){
    npart_x[ii] = -9999. ; 
    npart_y[ii] =  9999. ; 
  }

  for(int ii = 0; ii < 500; ii++){
    npart_x_of_nucleus_a[ii] = -9999. ; 
    npart_y_of_nucleus_a[ii] =  9999. ; 
    npart_w_of_nucleus_a[ii] =  0. ; 
  }

  for(int ii = 0; ii < 500; ii++){
    npart_x_of_nucleus_b[ii] = -9999. ; 
    npart_y_of_nucleus_b[ii] =  9999. ; 
    npart_w_of_nucleus_b[ii] =  0. ; 
  }

  for(int ii = 0; ii < 10000; ii++){
    ncoll_x[ii] = -9999. ; 
    ncoll_y[ii] =  9999. ; 
  }

  double sigma_perp_sq_ = sqr(sigma_perp) ;
  double smearing_range_ ;
  double contributors_x ; 
  double contributors_y ; 
  double contributors_w ; 

  // collect the npart and ncoll sources info from MC Glauber class.
  mc->get_npart_source_positions(npart_x, npart_y); 
  mc->get_npart_source_positions_of_nucleus_a(npart_x_of_nucleus_a, npart_y_of_nucleus_a); 
  mc->get_npart_source_positions_of_nucleus_b(npart_x_of_nucleus_b, npart_y_of_nucleus_b); 
  mc->get_npart_weight_of_nucleus_a(npart_w_of_nucleus_a); 
  mc->get_npart_weight_of_nucleus_b(npart_w_of_nucleus_b);   
  mc->get_ncoll_source_positions(ncoll_x, ncoll_y);
  // the sources are already adjusted to give CM = (0,0)
   
 
  // set the contribution to each cell to be zero at the beginning of the event.
  for(int ix = 0 ; ix < arena->get_nx() ; ix++ ){
    for(int iy = 0 ; iy < arena->get_ny() ; iy++ ){
      arena->get_cell(ix,iy)->set_contri_from_nucleus_a_after_gaussian_smearing(0.);
      arena->get_cell(ix,iy)->set_contri_from_nucleus_b_after_gaussian_smearing(0.);
      arena->get_cell(ix,iy)->set_contri_from_binary_collisions_after_gaussian_smearing(0.);
    }
  }


  // smear the participant sources of Nucleus A.
  for(int ipart=0; ipart < mc->get_no_of_participants_in_nucleus_a() ; ipart++ ){
    contributors_x = npart_x_of_nucleus_a[ipart] ; 
    contributors_y = npart_y_of_nucleus_a[ipart] ;
    contributors_w = npart_w_of_nucleus_a[ipart] ;
    smearing_range_ = 4. * sigma_perp ; 
    int upper_index_x, upper_index_y ; 
    get_nearest_cell_index( contributors_x + smearing_range_ , contributors_y + smearing_range_ , upper_index_x, upper_index_y ) ; 
    int lower_index_x, lower_index_y ;  
    get_nearest_cell_index( contributors_x - smearing_range_ , contributors_y - smearing_range_ , lower_index_x, lower_index_y ) ;
    for(int ix = lower_index_x ; ix < upper_index_x ; ix++ ){
      for(int iy = lower_index_y ; iy < upper_index_y ; iy++ ){
        double grid_x = -arena->get_xmax() + ix * arena->get_dx() ; 
        double grid_y = -arena->get_ymax() + iy * arena->get_dy() ;
        double prefactor_ = 1. / ( 2. * M_PI * sigma_perp_sq_ ) ; 
        double distance_sq = sqr(grid_x - contributors_x) + sqr(grid_y - contributors_y);
        double weight_a = prefactor_ * contributors_w * exp(-.5*distance_sq/sigma_perp_sq_); 
        arena->get_cell(ix,iy)->update_contri_from_nucleus_a_after_gaussian_smearing(weight_a);
      } // iy loop
    } // ix loop
  } // ipart loop


  // smear the participant sources of Nucleus B.
  for(int ipart=0; ipart < mc->get_no_of_participants_in_nucleus_b() ; ipart++ ){
    contributors_x = npart_x_of_nucleus_b[ipart] ; 
    contributors_y = npart_y_of_nucleus_b[ipart] ;
    contributors_w = npart_w_of_nucleus_b[ipart] ;
    smearing_range_ = 4. * sigma_perp ; 
    int upper_index_x, upper_index_y ; 
    get_nearest_cell_index( contributors_x + smearing_range_ , contributors_y + smearing_range_ , upper_index_x, upper_index_y ) ; 
    int lower_index_x, lower_index_y ;  
    get_nearest_cell_index( contributors_x - smearing_range_ , contributors_y - smearing_range_ , lower_index_x, lower_index_y ) ;
    for(int ix = lower_index_x ; ix < upper_index_x ; ix++ ){
      for(int iy = lower_index_y ; iy < upper_index_y ; iy++ ){
        double grid_x = -arena->get_xmax() + ix * arena->get_dx() ; 
        double grid_y = -arena->get_ymax() + iy * arena->get_dy() ;
        double prefactor_ = 1. / ( 2. * M_PI * sigma_perp_sq_ ) ; 
        double distance_sq = sqr(grid_x - contributors_x) + sqr(grid_y - contributors_y);
        double weight_b = prefactor_ * contributors_w * exp(-.5*distance_sq/sigma_perp_sq_); 
        arena->get_cell(ix,iy)->update_contri_from_nucleus_b_after_gaussian_smearing(weight_b);
      } // iy loop
    } // ix loop
  } // ipart loop


  /*
  // smear the binary collison sources.
  for(int icoll=0; icoll < mc->get_ncoll() ; icoll++ ){
    double contributors_x ; 
    double contributors_y ; 
    contributors_x = ncoll_x[icoll] ; 
    contributors_y = ncoll_y[icoll] ;
    int upper_index_x, upper_index_y ; 
    get_nearest_cell_index( contributors_x + 4 * sigma_perp , contributors_y + 4 * sigma_perp , upper_index_x, upper_index_y ) ; 
    int lower_index_x, lower_index_y ;  
    get_nearest_cell_index( contributors_x - 4 * sigma_perp , contributors_y - 4 * sigma_perp , lower_index_x, lower_index_y ) ;
    for(int ix = lower_index_x ; ix < upper_index_x ; ix++ ){
      for(int iy = lower_index_y ; iy < upper_index_y ; iy++ ){
        double grid_x = -arena->get_xmax() + ix * arena->get_dx() ; 
        double grid_y = -arena->get_ymax() + iy * arena->get_dy() ;
	   double weight_ncoll = 0 ; 
        arena->get_cell(ix,iy)->update_contri_from_binary_collisions_after_gaussian_smearing(weight_ncoll);
      } // iy loop
    } // ix loop
  } // icoll loop
  */
}


// relevant to generate Event-by-Event MC Glauber IC for hydro input. //
void mc_glau_smear::reset_contribution_to_zero_on_the_cells(){
  for(int ix = 0 ; ix < arena->get_nx() ; ix++ ){
    for(int iy = 0 ; iy < arena->get_ny() ; iy++ ){
      arena->get_cell(ix,iy)->reset_contributions_to_zero();
    }
  }
}


// Write the event averaged profile in MUSIC format.
void mc_glau_smear::write_event_averaged_profile_to_file_after_gaussian_smearing(int nEvents, int flag_to_generate_music_boost_invariant_file, int event_index){
  
  std::ofstream out_file;
  if(flag_to_generate_music_boost_invariant_file > 0 ){
    std::stringstream output_filename;
    output_filename.str("");      output_filename << "output/mc_glauber_single_event_transverse_profile_for_boost_invariant_music_" << event_index ;
    output_filename << ".dat";
    out_file.open(output_filename.str().c_str(), std::ios::out);
  }
  else{
    std::stringstream output_filename;
    output_filename.str("");
    output_filename << "output/mc_glauber_single_event_transverse_profile_for_rapidity_extension_" << event_index ;
    output_filename << ".dat";
    out_file.open(output_filename.str().c_str(), std::ios::out);
  }

  out_file <<"#"<<"\t"<<"mc_glauber"<<"\t"<<"1"<<"\t"<<"neta="<<"\t"<<"1"<<"\t"<<"nx="
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
        npart_contri_of_a = arena->get_cell(ix,iy)->get_contri_from_nucleus_a_after_gaussian_smearing() ; 
        npart_contri_of_b = arena->get_cell(ix,iy)->get_contri_from_nucleus_b_after_gaussian_smearing() ; 
        ncoll_contri      = arena->get_cell(ix,iy)->get_contri_from_binary_collisions_after_gaussian_smearing() ; 
      if((npart_contri_of_a+npart_contri_of_b) < 0.0000001 ){
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



double mc_glau_smear::get_reduced_thickness_after_gaussian_smearing() {
  double sum = 0.;
  double TA ;
  double TB ; 
  for (int iy = 0; iy < arena->get_ny(); ++iy) {
    for (int ix = 0; ix < arena->get_nx(); ++ix) {
      TA = arena->get_cell(ix,iy)->get_contri_from_nucleus_a_after_gaussian_smearing() ; 
      TB = arena->get_cell(ix,iy)->get_contri_from_nucleus_b_after_gaussian_smearing() ; 
      sum += reduced_thickness_func(TA,TB);
    }
  }
  return  arena->get_dx() * arena->get_dy() * sum;
}



// this function is taken from trento and modified
void mc_glau_smear::get_eccentricities_after_gaussian_smearing(double* eccentricity_) {
  // Compute eccentricity.

 // Simple helper class for use in the following loop.
  struct EccentricityAccumulator {
    double re = 0.;  // real part
    double im = 0.;  // imaginary part
    double wt = 0.;  // weight
    double finish() const  // compute final eccentricity
    { return std::sqrt(re*re + im*im) / std::fmax(wt, 1e-10); }
  } e2, e3, e4, e5;

  // first calculate xcm and ycm of the dist
  double sumx = 0.;
  double sumy = 0.;
  double sumw = 0.;
  double TA ;
  double TB ; 
  double TR ;
  double t ; 
  double cmx;
  double cmy;
  for (int iy = 0; iy < arena->get_ny(); ++iy) {
    for (int ix = 0; ix < arena->get_nx(); ++ix) {
      auto x = -arena->get_xmax() + ix * arena->get_dx() ;
      auto y = -arena->get_ymax() + iy * arena->get_dy() ;
      TA = arena->get_cell(ix,iy)->get_contri_from_nucleus_a_after_gaussian_smearing() ; 
      TB = arena->get_cell(ix,iy)->get_contri_from_nucleus_b_after_gaussian_smearing() ; 
      TR = reduced_thickness_func(TA,TB);
      sumx += x * TR ; 
      sumy += y * TR ; 
      sumw += TR ; 
    }
  }
  cmx = sumx / sumw ; 
  cmy = sumy / sumw ; 
 
 
  for (int iy = 0; iy < arena->get_ny(); ++iy) {
    for (int ix = 0; ix < arena->get_nx(); ++ix) {
      TA = arena->get_cell(ix,iy)->get_contri_from_nucleus_a_after_gaussian_smearing() ; 
      TB = arena->get_cell(ix,iy)->get_contri_from_nucleus_b_after_gaussian_smearing() ; 
      t = reduced_thickness_func(TA,TB);

      if (t < TINY)
        continue;

      // Compute (x, y) relative to the CM and cache powers of x, y, r.
      auto x = -arena->get_xmax() + ix * arena->get_dx() - cmx;
      auto x2 = x*x;
      auto x3 = x2*x;
      auto x4 = x2*x2;

      auto y = -arena->get_ymax() + iy * arena->get_dy() - cmy;
      auto y2 = y*y;
      auto y3 = y2*y;
      auto y4 = y2*y2;

      auto r2 = x2 + y2;
      auto r = std::sqrt(r2);
      auto r4 = r2*r2;

      auto xy = x*y;
      auto x2y2 = x2*y2;

      // The eccentricity harmonics are weighted averages of r^n*exp(i*n*phi)
      // over the entropy profile (reduced thickness).  The naive way to compute
      // exp(i*n*phi) at a given (x, y) point is essentially:
      //
      //   phi = arctan2(y, x)
      //   real = cos(n*phi)
      //   imag = sin(n*phi)
      //
      // However this implementation uses three unnecessary trig functions; a
      // much faster method is to express the cos and sin directly in terms of x
      // and y.  For example, it is trivial to show (by drawing a triangle and
      // using rudimentary trig) that
      //
      //   cos(arctan2(y, x)) = x/r = x/sqrt(x^2 + y^2)
      //   sin(arctan2(y, x)) = y/r = x/sqrt(x^2 + y^2)
      //
      // This is easily generalized to cos and sin of (n*phi) by invoking the
      // multiple angle formula, e.g. sin(2x) = 2sin(x)cos(x), and hence
      //
      //   sin(2*arctan2(y, x)) = 2*sin(arctan2(y, x))*cos(arctan2(y, x))
      //                        = 2*x*y / r^2
      //
      // Which not only eliminates the trig functions, but also naturally
      // cancels the r^2 weight.  This cancellation occurs for all n.
      //
      // The Event unit test verifies that the two methods agree.
      e2.re += t * (y2 - x2);
      e2.im += t * 2.*xy;
      e2.wt += t * r2;

      e3.re += t * (y3 - 3.*y*x2);
      e3.im += t * (3.*x*y2 - x3);
      e3.wt += t * r2*r;

      e4.re += t * (x4 + y4 - 6.*x2y2);
      e4.im += t * 4.*xy*(y2 - x2);
      e4.wt += t * r4;

      e5.re += t * y*(5.*x4 - 10.*x2y2 + y4);
      e5.im += t * x*(x4 - 10.*x2y2 + 5.*y4);
      e5.wt += t * r4*r;
    }
  }

  eccentricity_[2] = e2.finish();
  eccentricity_[3] = e3.finish();
  eccentricity_[4] = e4.finish();
  eccentricity_[5] = e5.finish();
}















