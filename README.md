   (GICG) Glauber-type Initial Condition Generator ver 3
   ====================================================
   
   Tribhuban Parida, PostDoc, AGH UST Krakow, Poland
   tparida@agh.edu.pl
   strting date : 27 Nov 2025(Thursday)  version 3
   This version is made to probe the v2-v3 fluctuation
   
   
   The code generates 2D/3D, smooth/fluctuating 
   initial profile using MC Glauber model which 
   can be used in hydro codes for evolution.

   eccentricity is calculated.

   // how to simulate
   1. Set input parameters in file "input.md".
   2. [command in terminal] make clean
   3. [command in terminal] make 
   4. [command in terminal] ./mc_glauber <no. of events want to generate>  <input.md>
   

  
 CHECKS TODO :
   1. check b distribution (minimum bias).  (Taking same Pb and U nucleus parameter as in trento)
   2. check Npart distribution (minimum bias). (Taking same Pb and U nucleus parameter as in trento)
   3. introduce gamma fluctuation and check whether gamma function is correct or not.    (checked )
   4. after introducing gamma fluctuation check the multiplicity distribution and compare with trento.
   5. Then compare e2, e3 distribution in different centrality and compare with trento.
   
   6. remove two component Glauber case.
   7. Recenter the energy distribution or not ? think about it !!!
