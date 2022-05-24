#  :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#  :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#  :::                                                                                           :::
#  :::  Set your input parameters here for MC-Glauber simulation.                                :::
#  :::  The code will not read any line of this file starting with symbol "#" or any empty line. :::
#  :::                                                                                           :::
#  :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#  :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

mode 1

#collision energy
#[Info] in GeV
SNN 200.0


#collision species & impact parameter
#[Info] species can be (Au,U,Pb,p) 
projectile Au2
target Au2
bmin  6.0001
bmax  6.0002

# Two component MC Glauber parametrs
npp 1.0
xhard 0.14

#minimum bias run root output file name, required for centrality determination.
root_output_file_name    output/Au_Au_200_min_bias_idnni_nx_201.root


#Grid, smearing and event averaged profile parameters
xmax 12
ymax 12
nx   241
ny   241
gaussian_smearing_sigma 0.4
upper_mult_proxy_cut  10000000.00
lower_mult_proxy_cut  0.000


#rapidity extension
etamax 8
neta   81
matter_eta_plateau 1.3
matter_eta_fall  1.5
matter_tilt_param  2.5
SARJ_model_f_param 0.25
init_time_tau0 0.8 

baryon_rapidity_profile_type 2

two_component_baryon_deposition_parameter_omega  0.1

baryon_rapidity_profile_eta_peak  4.0
baryon_rapidity_profile_sigma_eta_plus  0.6
baryon_rapidity_profile_sigma_eta_minus  0.8

iiserbpr_ra_type_2_peak  4.0
iiserbpr_ra_type_2_etal  1.0
iiserbpr_ra_type_2_right_fall  0.6



# :: END :: #






