#pragma once
#include<iostream>
#include<cmath>
#include <random>
#include <memory>
#define M_PI_ 3.1415927

class random_gen{

   private :
     long seed;
     std::random_device ran_dev;
     std::unique_ptr<std::mt19937> ran_generator;
     std::uniform_real_distribution<double> rand_uniform_dist;
     std::gamma_distribution<double> rand_gamma_dist;

   public :
     random_gen();
     random_gen(double k);
     ~random_gen();
     double rand_uniform();
     double rand_gamma();
     int get_seed();
     double random_gaussian(double, double);
     double random_sinx_0_to_pi();
     void random_3D_gaussian(double sigma, double pos[3]);
     double random_linear(double a, double b);
};
