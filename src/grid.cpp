#include "grid.h"

grid::grid(InputData *inparams_)
{
  inparams = inparams_ ; 
  
  // grid information below
  x_max = inparams->xmax ;
  y_max = inparams->ymax ;
  
  nx = inparams->nx ;
  ny = inparams->ny ;
  
  dx = ( 2 * x_max ) / ( nx - 1 ) ;
  dy = ( 2 * y_max ) / ( ny - 1 ) ;
  
  std::cout << "grid info below ... " << std::endl ; 
  std::cout << "x-max : " << x_max << "  y-max : " << y_max  << std::endl ;
  std::cout << "nx : " << nx << "  ny : " << ny  << std::endl ;
  std::cout << "dx : " << dx << "  dy : " << dy  << std::endl ;
  std::cout << "... ... ... ... ... " << std::endl ; 
  
  // create the grid with above properties
  create_grid();

}



grid::~grid()
{
  delete[] c ; 
}



void grid::create_grid()
{
  c = new cell[nx * ny];
  for (int ix = 0; ix < nx; ix++){
    for (int iy = 0; iy < ny; iy++){
      get_cell(ix, iy)->set_position(ix, iy);
    }
  }
}

