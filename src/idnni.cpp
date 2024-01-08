#include "idnni.h"

idnni::idnni(InputData *inparams_,grid* arena_, mc_glau* _mc)
{
  inparams = inparams_ ; 
  mc = _mc ;
  arena = arena_ ; 
  double radius_of_nucleon = sqrt ( mc->get_inelastic_cross_section() / (4*PI) )  ;
  
  // grid information below
  x_max = arena->get_xmax() ;
  y_max = arena->get_ymax() ;
  
  nx = arena->get_nx() ;
  ny = arena->get_ny() ;
  
  dx = arena->get_dx() ;
  dy = arena->get_dy() ;
  
  std::cout << "IDNNI grid info below ... " << std::endl ; 
  std::cout << "x-max : " << x_max << "  y-max : " << y_max  << std::endl ;
  std::cout << "nx : " << nx << "  ny : " << ny  << std::endl ;
  std::cout << "dx : " << dx << "  dy : " << dy  << std::endl ;
  std::cout << "radius of nucleons : " <<  radius_of_nucleon << std::endl ; 
  std::cout << "... ... ... ... ... " << std::endl ; 
  if( dx > radius_of_nucleon / 5  ){
    std::cout << "radius of nucleons : " <<  radius_of_nucleon << std::endl ; 
    std::cout << "dx > radius_of_nucleon / 5" << endl ; 
    std::cout << "make a fine grid ..." << std::endl ; 
    exit(1);
  }  
  
}


idnni::~idnni()
{
}


// Calculates the contribution on cells due to all nucleons in one event.
void idnni::perform()
{
  int A = mc->get_mass_number_of_nucleus_A();
  int B = mc->get_mass_number_of_nucleus_B();
  double radius_of_nucleon = sqrt ( mc->get_inelastic_cross_section() / (4*PI) )  ;
  
  
  mc->get_nucleus_A(XA, YA, ZA, pflagA);
  mc->get_nucleus_B(XB, YB, ZB, pflagB);
  
  // std::ofstream file_0;
  // file_0.open("nucleus_a.dat");
  // for(int i=0; i<A; i++)
  //    {
  //      file_0 << XA[i] << "\t" << YA[i] << "\t" << ZA[i] << endl; 
  //    }

  // file_0.close();
  // file_0.open("nucleus_b.dat");

  // for(int i=0; i<B; i++)
  //    {
  //      file_0 << XB[i] << "\t" << YB[i] << "\t" << ZB[i] << endl; 
  //    }
  // file_0.close();
  
  
  // Central index should be at '(nx+1)/2'.
  // But the grid index starts from zero in the code.
  // So in the code, the central index is at " (nx+1)/2 - 1 " .
  
  // set contributions to zero before starting one event.
  for(int ix=0; ix<nx; ix++){
    for(int iy=0; iy<ny; iy++){
      arena->get_cell(ix,iy)->set_contri_from_nucleus_a(0.);
      arena->get_cell(ix,iy)->set_contri_from_nucleus_b(0.);
    }
  }
  
  
  int share_n = 0 ;
  const int max_share_holders = ( 2 * radius_of_nucleon / dx ) * ( 2 * radius_of_nucleon / dy )  ; 
  int share_ix[max_share_holders] ;
  int share_iy[max_share_holders] ;
  for(int ii = 0 ; ii < max_share_holders ; ii++){
    share_ix[ii] = -5000 ;
    share_iy[ii] = -5000 ;  
  }
  
  
  // loop over nucleons of Nuceleus-A 
  for(int inucl = 0; inucl < A ; inucl ++){
    double participant_x = XA[inucl] ; 
    double participant_y = YA[inucl] ;
    
    int upper_index_x, upper_index_y ; 
    get_nearest_cell_index( participant_x + 2 * radius_of_nucleon , participant_y + 2 * radius_of_nucleon , upper_index_x, upper_index_y ) ; 
    int lower_index_x, lower_index_y ;  
    get_nearest_cell_index( participant_x - 2 * radius_of_nucleon , participant_y - 2 * radius_of_nucleon , lower_index_x, lower_index_y ) ; 
    // std::cout << "part-x : " << participant_x << "  part-y : " << participant_y << std::endl ; 
    // std::cout << "lx : " << -x_max + lower_index_x * dx << "  ly : " <<  -y_max + lower_index_y * dy << std::endl ; 
    // std::cout << "ux : " << -x_max + upper_index_x * dx << "  uy : " <<  -y_max + upper_index_y * dy << std::endl ;
    
    share_n = 0 ; 
    for(int ix = lower_index_x ; ix < upper_index_x ; ix++ ){
      for(int iy = lower_index_y ; iy < upper_index_y ; iy++ ){
	double grid_x = -x_max + ix * dx ; 
	double grid_y = -y_max + iy * dy ;
	double distance = sqrt ( pow( participant_x-grid_x, 2)  + pow( participant_y-grid_y, 2 ) ) ;
	if (distance < radius_of_nucleon ){ // if one cell falls within the radius of nucleon then consider that as a share holder.
	  share_ix[share_n] = ix ;  // store the ix-index of the share holder.
	  share_iy[share_n] = iy ;  // store the iy-index of the share holder.
	  share_n += 1 ;            // incrase the number of share holders by one unit. 
	}
      } // iy loop
    } // ix loop
    
    for(int ii = 0 ; ii < share_n ; ii++ ){
      //std::cout << "         ......................       " << std::endl ;
      //std::cout << "share number : " << ii << std::endl ; 
      //std::cout << "participant @(" << participant_x << ", " << participant_y << ")" << std::endl ; 
      //std::cout << "cell (" << -x_max + share_ix[ii] * dx << ", " << -y_max + share_iy[ii] * dy << ")  share :" <<  1.0 / share_n  << std::endl ;  
      arena->get_cell(share_ix[ii],share_iy[ii])->update_contri_from_nucleus_a( 1.0 / share_n ) ;  // Distribute the shares among share holders equally. 
    }
    
  } // inucl loop
  
  
  
  // loop over nucleons of Nuceleus-B 
  for(int inucl = 0; inucl < B ; inucl ++){
    double participant_x = XB[inucl] ; 
    double participant_y = YB[inucl] ;
    
    int upper_index_x, upper_index_y ; 
    get_nearest_cell_index( participant_x + 2 * radius_of_nucleon , participant_y + 2 * radius_of_nucleon , upper_index_x, upper_index_y ) ; 
    int lower_index_x, lower_index_y ;  
    get_nearest_cell_index( participant_x - 2 * radius_of_nucleon , participant_y - 2 * radius_of_nucleon , lower_index_x, lower_index_y ) ; 
    
    share_n = 0 ; 
    for(int ix = lower_index_x ; ix < upper_index_x ; ix++ ){
      for(int iy = lower_index_y ; iy < upper_index_y ; iy++ ){
	double grid_x = -x_max + ix * dx ; 
	double grid_y = -y_max + iy * dy ;
	double distance = sqrt ( pow( participant_x-grid_x, 2)  + pow( participant_y-grid_y, 2 ) ) ;
	if (distance < radius_of_nucleon ){
	  share_ix[share_n] = ix ; 
	  share_iy[share_n] = iy ; 
	  share_n += 1 ; 
	}
      } // iy loop
    } // ix loop
    
    for(int ii = 0 ; ii < share_n ; ii++ ){
      arena->get_cell(share_ix[ii],share_iy[ii])->update_contri_from_nucleus_b( 1.0 / share_n ) ; 
    }
    
  } // inucl loop
  
}



void idnni::get_nearest_cell_index(double xx, double yy, int& ix, int& iy ){
  ix =  ceil( ( xx - (-x_max) ) / dx ) ; 
  iy =  ceil( ( yy - (-y_max) ) / dy ) ; 
}


// calculates the effective Npart.
double idnni::get_total_deposition_over_all_cells_for_current_event(){
  double net_dist = 0 ; 
  for(int ix = 0 ; ix < nx ; ix++ ){
    for(int iy = 0 ; iy < ny ; iy++ ){
      if( arena->get_cell(ix,iy)->get_contri_from_nucleus_a() > 0.000001 && arena->get_cell(ix,iy)->get_contri_from_nucleus_b() > 0.000001 )
	net_dist += arena->get_cell(ix,iy)->get_contri_from_nucleus_a() + arena->get_cell(ix,iy)->get_contri_from_nucleus_b() ; 
    }
  }
  return net_dist ; 
}


void idnni::get_eccentricities_and_participant_plane_angles(int Norder, double* eps, double* psi){
  for(int iorder = 0 ; iorder < Norder ; iorder++ ){
    eps[iorder] = 0. ; 
    psi[iorder] = 0. ; 
  }
  
  double cosnphi[Norder] ;  // represents <cos(n*phi)>
  double sinnphi[Norder] ;  // represents <sin(n*phi)>
  double rn_dot_depo[Norder] ; 
  for(int iorder =0 ; iorder < Norder ; iorder++ ){
    cosnphi[iorder] = 0. ; 
    sinnphi[iorder] = 0. ; 
    rn_dot_depo[iorder] = 0. ; 
  }
  
  // (xx,yy) is the center of mass of the distribution.
  double xx,yy ; 
  get_cm_of_the_distribution(xx,yy);
  
  for(int ix = 0 ; ix < nx ; ix++ ){
    for(int iy = 0 ; iy < ny ; iy++ ){
      if( arena->get_cell(ix,iy)->get_contri_from_nucleus_a() < 0.000001 || arena->get_cell(ix,iy)->get_contri_from_nucleus_b() < 0.000001 )
	continue ;
      double grid_x = -x_max + ix * dx ; 
      double grid_y = -y_max + iy * dy ;
      grid_x -= xx ; // shifting the origin to CM
      grid_y -= yy ; 
      double r = sqrt(grid_x * grid_x + grid_y * grid_y) ; 
      double phi = atan2(grid_y, grid_x) ; 
      for(int iorder =0 ; iorder < Norder ; iorder++ ){
	if(iorder == 1 ){
	  cosnphi[iorder] += pow(r, 3) * cos(iorder*phi) * arena->get_cell(ix,iy)->get_net_deposition_on_the_cell() ; 
	  sinnphi[iorder] += pow(r, 3) * sin(iorder*phi) * arena->get_cell(ix,iy)->get_net_deposition_on_the_cell() ;
	  rn_dot_depo[iorder] += pow(r, 3) * arena->get_cell(ix,iy)->get_net_deposition_on_the_cell() ; // sum over ( r^{n} * weight )
	}
	else{
	  cosnphi[iorder] += pow(r, iorder) * cos(iorder*phi) * arena->get_cell(ix,iy)->get_net_deposition_on_the_cell() ; 
	  sinnphi[iorder] += pow(r, iorder) * sin(iorder*phi) * arena->get_cell(ix,iy)->get_net_deposition_on_the_cell() ; 
          rn_dot_depo[iorder] += pow(r, iorder) * arena->get_cell(ix,iy)->get_net_deposition_on_the_cell() ; 
	}
      } // iorder loop
    } // iy loop
  } // ix loop
  
  // <cos(0*phi)> is 1.
  // <sin(0*phi)> is 0.
  for(int iorder = 1 ; iorder < Norder ; iorder++ ){
    cosnphi[iorder] /= -rn_dot_depo[iorder] ; 
    sinnphi[iorder] /= -rn_dot_depo[iorder] ; 
  }
  
  for(int iorder = 1 ; iorder < Norder ; iorder++ ){
    eps[iorder] = sqrt( cosnphi[iorder] * cosnphi[iorder] + sinnphi[iorder] * sinnphi[iorder] )  ; 
    psi[iorder] = atan2(sinnphi[iorder],  cosnphi[iorder]) /  iorder ; 
  }
}


void idnni::get_cm_of_the_distribution(double& xx, double& yy){
  xx = -10000.0 ; 
  yy = -10000.0 ;
  double density_sum = 0. ;
  double x_dot_density_sum = 0. ;  
  double y_dot_density_sum = 0. ;  
  for(int ix = 0 ; ix < nx ; ix++ ){
    for(int iy = 0 ; iy < ny ; iy++ ){
      if( arena->get_cell(ix,iy)->get_contri_from_nucleus_a() < 0.000001 || arena->get_cell(ix,iy)->get_contri_from_nucleus_b() < 0.000001 )
	continue ; 
      double grid_x = -x_max + ix * dx ; 
      double grid_y = -y_max + iy * dy ;
      density_sum += arena->get_cell(ix,iy)->get_net_deposition_on_the_cell() ;
      x_dot_density_sum += grid_x * arena->get_cell(ix,iy)->get_net_deposition_on_the_cell();
      y_dot_density_sum += grid_y * arena->get_cell(ix,iy)->get_net_deposition_on_the_cell();
    }
  }
  xx = x_dot_density_sum / density_sum ;  
  yy = y_dot_density_sum / density_sum ;  
}



// * * * * * * * * * * * * * * * * * * * * * *
// Below part is for event averaged profile  *
// * * * * * * * * * * * * * * * * * * * * * *


// Event by Event the contributions are added to the cells. 
// Later averaging on cells over all events will take place.
void idnni::update_contribution_on_cells_over_all_events(){
  for(int ix = 0 ; ix < nx ; ix++ ){
    for(int iy = 0 ; iy < ny ; iy++ ){
      if( arena->get_cell(ix,iy)->get_contri_from_nucleus_a() > 0.000001 && arena->get_cell(ix,iy)->get_contri_from_nucleus_b() > 0.000001 ){
	arena->get_cell(ix,iy)->update_contribution_from_nucleus_a_over_all_events(arena->get_cell(ix,iy)->get_contri_from_nucleus_a()) ;
	arena->get_cell(ix,iy)->update_contribution_from_nucleus_b_over_all_events(arena->get_cell(ix,iy)->get_contri_from_nucleus_b()) ;
      }
    }
  }
}



// * * * * * * * * * * * * * * * * * * * * * *
//   Below part is for gaussian smearing     *
// * * * * * * * * * * * * * * * * * * * * * *

void idnni::update_contribution_on_cells_over_all_events_with_gaussian_smearing(){
  for(int ix = 0 ; ix < nx ; ix++ ){
    for(int iy = 0 ; iy < ny ; iy++ ){
      if( arena->get_cell(ix,iy)->get_contri_from_nucleus_a_after_gaussian_smearing() > 0.000001 && 
	  arena->get_cell(ix,iy)->get_contri_from_nucleus_b_after_gaussian_smearing() > 0.000001 ){
	arena->get_cell(ix,iy)->update_contribution_from_nucleus_a_over_all_events(arena->get_cell(ix,iy)->get_contri_from_nucleus_a_after_gaussian_smearing()) ;
	arena->get_cell(ix,iy)->update_contribution_from_nucleus_b_over_all_events(arena->get_cell(ix,iy)->get_contri_from_nucleus_b_after_gaussian_smearing()) ;
      }
    }
  }
}




void idnni::do_gaussian_smearing_for_the_current_event(double sigma_perp){
  // double variance = sigma_perp * sigma_perp ; 
  for(int ix = 0 ; ix < nx ; ix++ ){
    for(int iy = 0 ; iy < ny ; iy++ ){
      arena->get_cell(ix,iy)->set_contri_from_nucleus_a_after_gaussian_smearing(0.);
      arena->get_cell(ix,iy)->set_contri_from_nucleus_b_after_gaussian_smearing(0.);
    }
  }
  
  for(int ix = 0 ; ix < nx ; ix++ ){
    for(int iy = 0 ; iy < ny ; iy++ ){
      // going to findout the contribution by other cells to the current cell(ix,iy).
      double grid_x = -x_max + ix * dx ; 
      double grid_y = -y_max + iy * dy ;
      int upper_index_x, upper_index_y ; 
      get_nearest_cell_index( grid_x + 4 * sigma_perp , grid_y + 4 * sigma_perp , upper_index_x, upper_index_y ) ; 
      int lower_index_x, lower_index_y ;  
      get_nearest_cell_index( grid_x - 4 * sigma_perp , grid_y - 4 * sigma_perp , lower_index_x, lower_index_y ) ; 
      for(int jx = lower_index_x ; jx < upper_index_x ; jx++ ){
	for(int jy = lower_index_y ; jy < upper_index_y ; jy++ ){
	  if( arena->get_cell(jx,jy)->get_contri_from_nucleus_a() < 0.000001 || arena->get_cell(jx,jy)->get_contri_from_nucleus_b() < 0.000001  )
	    continue ; 
          double contributors_x ; 
          double contributors_y ; 
	  contributors_x = -x_max + jx * dx ; 
	  contributors_y = -y_max + jy * dy ;
	  double weight_a = 0.25 * ( -erf( (grid_x - dx / 2 - contributors_x) / sqrt(2) / sigma_perp ) + 
				     erf((grid_x + dx / 2 - contributors_x) / sqrt(2) / sigma_perp) )
	    * ( -erf( (grid_y - dy / 2 - contributors_y) / sqrt(2) / sigma_perp ) + 
		erf((grid_y + dy / 2 - contributors_y) / sqrt(2) / sigma_perp) )
	    * arena->get_cell(jx,jy)->get_contri_from_nucleus_a() ; 
	  //double weight_a = 1. / ( 2 * PI * variance ) *
	  //         exp( - pow(grid_x - contributors_x, 2) / (2 * variance) - pow(grid_y - contributors_y, 2) / (2 * variance) ) *
	  //            arena->get_cell(jx,jy)->get_contri_from_nucleus_a() * dx * dy ;
	  arena->get_cell(ix,iy)->update_contri_from_nucleus_a_after_gaussian_smearing(weight_a);
	  double weight_b = 0.25 * ( -erf( (grid_x - dx / 2 - contributors_x) / sqrt(2) / sigma_perp ) + 
				     erf((grid_x + dx / 2 - contributors_x) / sqrt(2) / sigma_perp) )
	    * ( -erf( (grid_y - dy / 2 - contributors_y) / sqrt(2) / sigma_perp ) + 
		erf((grid_y + dy / 2 - contributors_y) / sqrt(2) / sigma_perp) )
	    * arena->get_cell(jx,jy)->get_contri_from_nucleus_b() ; 
	  //double weight_b = 1. / ( 2 * PI * variance ) *
	  //         exp( - pow(grid_x - contributors_x, 2) / (2 * variance) - pow(grid_y - contributors_y, 2) / (2 * variance) ) *
	  //            arena->get_cell(jx,jy)->get_contri_from_nucleus_b() * dx * dy ;
	  arena->get_cell(ix,iy)->update_contri_from_nucleus_b_after_gaussian_smearing(weight_b);
	}
      }
    }
  }
} 




void idnni::rotation_by_second_order_event_plane_followed_by_gaussian_smearing_for_the_current_event(double sigma_perp, double rotation_angle){
  // double variance = sigma_perp * sigma_perp ; 
  for(int ix = 0 ; ix < nx ; ix++ ){
    for(int iy = 0 ; iy < ny ; iy++ ){
      arena->get_cell(ix,iy)->set_contri_from_nucleus_a_after_gaussian_smearing(0.);
      arena->get_cell(ix,iy)->set_contri_from_nucleus_b_after_gaussian_smearing(0.);
    }
  }

  for(int ix = 0 ; ix < nx ; ix++ ){
    for(int iy = 0 ; iy < ny ; iy++ ){
      double grid_x = -x_max + ix * dx ; 
      double grid_y = -y_max + iy * dy ;
      if( arena->get_cell(ix,iy)->get_contri_from_nucleus_a() < 0.000001 || arena->get_cell(ix,iy)->get_contri_from_nucleus_b() < 0.000001  )
        continue ; 
      double contributors_x ; // a grid point is a contributor now
      double contributors_y ; 
      contributors_x =  cos(rotation_angle) * grid_x + sin(rotation_angle) * grid_y ; 
      contributors_y = -sin(rotation_angle) * grid_x + cos(rotation_angle) * grid_y ;
      int upper_index_x, upper_index_y ; 
      get_nearest_cell_index( contributors_x + 4 * sigma_perp , contributors_y + 4 * sigma_perp , upper_index_x, upper_index_y ) ; 
      int lower_index_x, lower_index_y ;  
      get_nearest_cell_index( contributors_x - 4 * sigma_perp , contributors_y - 4 * sigma_perp , lower_index_x, lower_index_y ) ;
      if(upper_index_x > nx ){cout << "out of grid after rotation." << endl; exit(1); } 
      if(upper_index_y > ny ){cout << "out of grid after rotation." << endl; exit(1); } 
      if(lower_index_x < 0 ) {cout << "out of grid after rotation." << endl; exit(1); } 
      if(lower_index_y < 0 ) {cout << "out of grid after rotation." << endl; exit(1); } 

      for(int jx = lower_index_x ; jx < upper_index_x ; jx++ ){
	for(int jy = lower_index_y ; jy < upper_index_y ; jy++ ){
      double local_x = -x_max + jx * dx ; 
      double local_y = -y_max + jy * dy ;
	  double weight_a = 0.25 * ( -erf( (local_x - dx / 2 - contributors_x) / sqrt(2) / sigma_perp ) + 
				     erf((local_x + dx / 2 - contributors_x) / sqrt(2) / sigma_perp) )
	    * ( -erf( (local_y - dy / 2 - contributors_y) / sqrt(2) / sigma_perp ) + 
		erf((local_y + dy / 2 - contributors_y) / sqrt(2) / sigma_perp) )
	    * arena->get_cell(ix,iy)->get_contri_from_nucleus_a() ;
	  arena->get_cell(jx,jy)->update_contri_from_nucleus_a_after_gaussian_smearing(weight_a); 
          
	  double weight_b = 0.25 * ( -erf( (local_x - dx / 2 - contributors_x) / sqrt(2) / sigma_perp ) + 
				     erf((local_x + dx / 2 - contributors_x) / sqrt(2) / sigma_perp) )
	    * ( -erf( (local_y - dy / 2 - contributors_y) / sqrt(2) / sigma_perp ) + 
		erf((local_y + dy / 2 - contributors_y) / sqrt(2) / sigma_perp) )
	    * arena->get_cell(ix,iy)->get_contri_from_nucleus_b() ;
	  arena->get_cell(jx,jy)->update_contri_from_nucleus_b_after_gaussian_smearing(weight_b);      
        } // jy loop
      } // jx loop
    }
  }
}



void idnni::get_cm_of_the_distribution_after_gaussian_smearing(double& xx, double& yy){
  xx = -10000.0 ; 
  yy = -10000.0 ;
  double density_sum = 0. ;
  double x_dot_density_sum = 0. ;  
  double y_dot_density_sum = 0. ;  
  for(int ix = 0 ; ix < nx ; ix++ ){
    for(int iy = 0 ; iy < ny ; iy++ ){
      if( arena->get_cell(ix,iy)->get_contri_from_nucleus_a_after_gaussian_smearing() < 0.000001 || 
          arena->get_cell(ix,iy)->get_contri_from_nucleus_b_after_gaussian_smearing() < 0.000001 ){
	continue ; 
      }
      double grid_x = -x_max + ix * dx ; 
      double grid_y = -y_max + iy * dy ;
      density_sum += arena->get_cell(ix,iy)->get_net_deposition_on_the_cell_after_gaussian_smearing() ;
      x_dot_density_sum += grid_x * arena->get_cell(ix,iy)->get_net_deposition_on_the_cell_after_gaussian_smearing();
      y_dot_density_sum += grid_y * arena->get_cell(ix,iy)->get_net_deposition_on_the_cell_after_gaussian_smearing();
    }
  }
  xx = x_dot_density_sum / density_sum ;  
  yy = y_dot_density_sum / density_sum ;  
}


// calculates the effective Npart after gaussian smearing.
double idnni::get_total_deposition_over_all_cells_after_gaussian_smearing_for_current_event(){
  double net_dist = 0 ; 
  for(int ix = 0 ; ix < nx ; ix++ ){
    for(int iy = 0 ; iy < ny ; iy++ ){
      if( arena->get_cell(ix,iy)->get_contri_from_nucleus_a_after_gaussian_smearing() > 0.000001 && 
	  arena->get_cell(ix,iy)->get_contri_from_nucleus_b_after_gaussian_smearing() > 0.000001 ){
	net_dist += arena->get_cell(ix,iy)->get_contri_from_nucleus_a_after_gaussian_smearing() + 
	  arena->get_cell(ix,iy)->get_contri_from_nucleus_b_after_gaussian_smearing() ; 
      }
    }
  }
  return net_dist ; 
}


void idnni::get_eccentricities_and_participant_plane_angles_after_gaussian_smearing(int Norder, double* eps, double* psi){
  for(int iorder = 0 ; iorder < Norder ; iorder++ ){
    eps[iorder]  = 0. ; 
    psi[iorder]  = 0. ; 
  }
  
  double cosnphi[Norder] ;  // represents <cos(n*phi)>
  double sinnphi[Norder] ;  // represents <sin(n*phi)>
  double rn_dot_depo[Norder] ; 
  for(int iorder =0 ; iorder < Norder ; iorder++ ){
    cosnphi[iorder] = 0. ; 
    sinnphi[iorder] = 0. ; 
    rn_dot_depo[iorder] = 0. ; 
  }
  
  // (xx,yy) is the center of mass of the distribution.
  double xx,yy ; 
  get_cm_of_the_distribution_after_gaussian_smearing(xx,yy);
  
  for(int ix = 0 ; ix < nx ; ix++ ){
    for(int iy = 0 ; iy < ny ; iy++ ){
      if( arena->get_cell(ix,iy)->get_contri_from_nucleus_a_after_gaussian_smearing() < 0.000001 || 
	  arena->get_cell(ix,iy)->get_contri_from_nucleus_b_after_gaussian_smearing() < 0.000001 ){
	continue ;
      }
      double grid_x = -x_max + ix * dx ; 
      double grid_y = -y_max + iy * dy ;
      grid_x -= xx ; // shifting the origin to CM
      grid_y -= yy ; 
      double r = sqrt(grid_x * grid_x + grid_y * grid_y) ; 
      double phi = atan2(grid_y, grid_x) ; 
      for(int iorder =0 ; iorder < Norder ; iorder++ ){
	if(iorder == 1 ){
	  cosnphi[iorder] += pow(r, 3) * cos(iorder*phi) * arena->get_cell(ix,iy)->get_net_deposition_on_the_cell_after_gaussian_smearing() ; 
	  sinnphi[iorder] += pow(r, 3) * sin(iorder*phi) * arena->get_cell(ix,iy)->get_net_deposition_on_the_cell_after_gaussian_smearing() ;
	  rn_dot_depo[iorder] += pow(r, 3) * arena->get_cell(ix,iy)->get_net_deposition_on_the_cell_after_gaussian_smearing() ; // sum over ( r^{n} * weight )
	}
	else{
	  cosnphi[iorder] += pow(r, iorder) * cos(iorder*phi) * arena->get_cell(ix,iy)->get_net_deposition_on_the_cell_after_gaussian_smearing() ; 
	  sinnphi[iorder] += pow(r, iorder) * sin(iorder*phi) * arena->get_cell(ix,iy)->get_net_deposition_on_the_cell_after_gaussian_smearing() ; 
          rn_dot_depo[iorder] += pow(r, iorder) * arena->get_cell(ix,iy)->get_net_deposition_on_the_cell_after_gaussian_smearing() ; 
	}
      } // iorder loop
    } // iy loop
  } // ix loop
  
  // <cos(0*phi)> is 1.
  // <sin(0*phi)> is 0.
  for(int iorder = 1 ; iorder < Norder ; iorder++ ){
    cosnphi[iorder] /= -rn_dot_depo[iorder] ; 
    sinnphi[iorder] /= -rn_dot_depo[iorder] ; 
  }
  
  for(int iorder = 1 ; iorder < Norder ; iorder++ ){
    eps[iorder] = sqrt( cosnphi[iorder] * cosnphi[iorder] + sinnphi[iorder] * sinnphi[iorder] )  ; 
    psi[iorder] = atan2(sinnphi[iorder],  cosnphi[iorder]) / iorder ; 
  }
}



// Write the event averaged profile in MUSIC format.
void idnni::write_event_averaged_profile_to_file_after_gaussian_smearing(int nEvents, int flag_to_generate_music_boost_invariant_file){
  
  std::ofstream out_file;
  if(flag_to_generate_music_boost_invariant_file > 0 ){
    out_file.open("output/idnni_event_averaged_profile_for_boost_invariant_music.dat", std::ios::out);
    cout << "the output file could be directly used in boost invariant music ... " << endl ; 
  }
  else{
    out_file.open("output/idnni_boost_invariant_event_averaged_profile_for_rapidity_extension.dat", std::ios::out);
  }

  out_file<<"#"<<"\t"<<"event_avergaed_glauber"<<"\t"<<"1"<<"\t"<<"neta="<<"\t"<<"1"<<"\t"<<"nx="<<"\t"<<nx<<"\t"<<"ny="<<"\t"<<ny
	  <<"\t"<<"deta="<<"\t"<<"0.1"<<"\t"<<"dx="<<"\t"<<dx<<"\t"<<"dy="<<"\t"<<dy<<endl;
  
  double mult   = 0.0 ;
  double mult_a = 0.0 ;
  double mult_b = 0.0 ;
  double participant_sum  = 0. ; 
  double participant_diff = 0. ; 

  for(int ix = 0 ; ix < nx ; ix++ ){
    for(int iy = 0 ; iy < ny ; iy++ ){
      double grid_x = -x_max + ix * dx ; 
      double grid_y = -y_max + iy * dy ;
      if( arena->get_cell(ix,iy)->get_contribution_from_nucleus_a_over_all_events() > 0.000001 && 
	  arena->get_cell(ix,iy)->get_contribution_from_nucleus_b_over_all_events() > 0.000001 ){
	mult = arena->get_cell(ix,iy)->get_contribution_from_nucleus_a_over_all_events() +
	  arena->get_cell(ix,iy)->get_contribution_from_nucleus_b_over_all_events() ; 
        mult_a = arena->get_cell(ix,iy)->get_contribution_from_nucleus_a_over_all_events() ; 
        mult_b = arena->get_cell(ix,iy)->get_contribution_from_nucleus_b_over_all_events() ; 
      }
      else{
	mult   = 0.000000002 ; 
	mult_a = 0.000000001 ; 
	mult_b = 0.000000001 ; 
      }

      mult    /=  nEvents ; 
      mult_a  /=  nEvents ; 
      mult_b  /=  nEvents ; 

      participant_sum += mult ;
      participant_diff += (mult_a - mult_b) ;  
      
      if(flag_to_generate_music_boost_invariant_file  > 0 ){ 
        out_file << "0" << "\t" << grid_x << "\t" << grid_y << "\t" << mult   
               << "\t" << "1" << "\t" << "0" << "\t" << "0" << "\t" << "0"
               << "\t" << "0" << "\t" << "0" << "\t" << "0" << endl ;
      }
      else{
        out_file << "0" << "\t" << grid_x << "\t" << grid_y << "\t" << mult << "\t" << mult_a  << "\t" << mult_b  
               << "\t" << "1" << "\t" << "0" << "\t" << "0" << "\t" << "0"
               << "\t" << "0" << "\t" << "0" << "\t" << "0" << endl ;
      }  
    } // iy loop
  } // ix loop

  // double ybeam = acosh(inparams->SNN / ( 2. * 0.938 ) ) ; 
  cout << "participant assymetry : " << ( participant_diff / participant_sum )  * 100   << endl ; 

}


