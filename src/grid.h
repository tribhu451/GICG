#pragma once
#include <fstream>
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <TRandom3.h>
#include <TF1.h>
#include "TMath.h"
#include "input_data.h"
#include "mc_glau.h"
#include "cell.h"
#include <gsl/gsl_sf_erf.h>


#define PI  3.1415927

using std::cout;
using std::endl;


class grid{

public : 

  grid(InputData *inparams_);
  ~grid();

  inline cell *get_cell(int ix, int iy){
    ix = ix > 0 ? ix : 0;
    ix = ix < nx ? ix : nx - 1;
    iy = iy > 0 ? iy : 0;
    iy = iy < ny ? iy : ny - 1;
    return &c[ix + nx * iy ];
  }
  
  void set_xmax(double xx){x_max = xx;}
  void set_ymax(double xx){y_max = xx;}
  void set_nx(int xx){nx = xx;}
  void set_ny(int xx){ny = xx;}
  void set_dx(double xx){dx = xx ;}
  void set_dy(double xx){dy = xx ;}


  inline double get_xmax(){ return x_max ;}
  inline double get_ymax(){ return y_max;}
  inline int    get_nx()   { return nx;}
  inline int    get_ny()   { return ny;}
  inline double get_dx()   { return dx;}
  inline double get_dy()   { return dy;}


private :
  cell* c; 
  InputData *inparams;
  void create_grid();
  // grid info
  float x_max;
  float y_max;
  int nx;
  int ny;
  float dx;
  float dy;
  void get_nearest_cell_index(double xx, double yy, int& ix, int& iy ) ; 

};
