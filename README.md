   =====================================================================
   ==   (GICG) Glauber-type Initial Condition Generator ver 0.0.0.2   ==
   =====================================================================

   Please inform authors before using this version ...
   Email : tribhu.451@gmail.com
      
   In this version EM filed is generated from Glauber spectator config.
   The current version could only work for
   1. event averged simulations
   2. For symmetric nucleus collision (simple modification is feasible for asymmetric ones.)

   Tribhuban Parida, Research Scholar, IISER Berhampur
   tribhubanp18@iiserbpr.ac.in
   strting date : Dec 2023  version 0.0.0.2
   
   The code generates smooth 
   initial profile using MC Glauber model which 
   can be used in hydro codes for evolution.

   eccentricity is calculated assuming two 
   component energy deposition with
   xhard in MC Glauber model .

   // how to simulate
   1. Set input parameters in file "input.md".
   2. [command in terminal] make clean
   3. [command in terminal] make 
   4. [command in terminal] ./mc_glauber_spc  <no. of events want to generate>  <input.md>
   

  

