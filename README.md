   (GICG) Glauber-type Initial Condition Generator ver 3
   ====================================================
   
   Tribhuban Parida, PostDoc, AGH UST Krakow, Poland
   Email : tparida@agh.edu.pl
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
   1. after introducing gamma fluctuation check the multiplicity distribution and compare with trento.  
     (Taking same Pb and U nucleus parameter as in trento)
   2. Then compare e2, e3 distribution in different centrality and compare with trento.

   The above checks are already done and you could find the corresponding plots at : 
   https://drive.google.com/drive/folders/1EaLLZDX34hJFsxch9jTUaxnjELUafcWZ


