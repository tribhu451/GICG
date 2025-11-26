#include "random.h"



random_gen::random_gen() : rand_uniform_dist(0.,1.) { 
  seed = ran_dev();
  //seed = 12345 ;  // avoid the randomness by turning on this one (only for check purpose) ... :)
  ran_generator = std::unique_ptr<std::mt19937>(new std::mt19937(seed));

}

random_gen::~random_gen(){
}


double random_gen::rand_uniform(){
  return(rand_uniform_dist(*ran_generator));
}


int random_gen::get_seed() {
  return(seed);
}


// This function generates the random numbers with gaussian 
// distribution having mean=mu and std. dev. = sigma
double random_gen::random_gaussian(double mu, double sigma){
  double x = 0. ; 

  if(mu<0.00000){
    sigma *= -1 ; 
  }
  double limit_of_call = mu + 4.0 * sigma ; 

  int ii=0 ; 
  for( int i=0; i<1000; i++ ){
    ii += 1 ; 
    double test_rand = rand_uniform();
    x = 2.0 * limit_of_call * rand_uniform() - limit_of_call ;
    double gauss = exp( -0.5 * pow( (x - mu) / sigma , 2 ) ) ;  
    if( test_rand < gauss ){
      break ; 
    }
  }

  if(ii > 997 ){
    std::cout << "maximum calls for normal" 
              << "distribution sampling ... " << std::endl ;
    std::cout << "exiting ... " << std::endl ;  
    exit(1);
  }
  
  return x ; 
 
}


// This function generates a random number
// between 0 to Pi with sin(x) probability
double random_gen::random_sinx_0_to_pi(){
  double cos_theta = 2 * rand_uniform() - 1  ; 
  double sin_theta = sqrt( 1.0 - cos_theta * cos_theta );
  int aa =1;
  double rrx = 2 * rand_uniform() - 1 ; 

  if(rrx<0){
    aa = 1 ; 
  }else{
    aa = -1 ; 
  }
   
  double thet = asin(aa*sin_theta)  ;
  if (thet<0) {thet += M_PI;}
  return  thet; 
}

// random number with P(x) ~ x in [a, b]
double random_gen::random_linear(double a, double b) {
    double r = rand_uniform();      // uniform random number [0,1]
    return a + (b - a) * sqrt(r);  // inverse CDF
}



// This function generates 3 random numbers (x,y,z) within a 3D Gaussian profile
// where the radius r = sqrt(x^2+y^2+z^2) follows a gaussian distribution 
// having mean=0 and std. dev. = sigma.
void random_gen::random_3D_gaussian(double sigma, double pos[3]){
  double r = 404 ; 
  double limit_of_call = 4.0 * sigma ; 

  int ii=0 ; 
  for( int i=0; i<10000; i++ ){
    ii += 1 ; 
    double test_rand = rand_uniform();
    r =  limit_of_call * rand_uniform()  ;
    double gauss = exp( -1.5 * pow( r / sigma , 2 ) ) ;  
    if( test_rand < gauss ){
      break ; 
    }
  }

  if(ii > 9997 ){
    std::cout << "maximum calls for normal" 
              << "distribution sampling ... " << std::endl ;
    std::cout << "exiting ... " << std::endl ;  
    exit(-1);
  }
   
  double th = random_sinx_0_to_pi();
  double phi = 2*M_PI*rand_uniform();

  pos[0] = r * sin(th) * cos(phi); // x
  pos[1] = r * sin(th) * sin(phi); // y
  pos[2] = r * cos(th); // z
 
}
