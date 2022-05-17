#  :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#  :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#  :::                                                                                           :::
#  :::  Set your input parameters here for MC-Glauber simulation.                                :::
#  :::  The code will not read any line of this file starting with symbol "#" or any empty line. :::
#  :::                                                                                           :::
#  :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#  :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::



#collision energy
#[Info] in GeV
SNN 200.0


#collision species & impact parameter
#[Info] species can be (Au,U,Pb,p) 
projectile Au2
target Au2
bmin  0.0001
bmax  10.00

npp 2.25
xhard 0.14

#output file name
root_output_file_name    output/Au_Au_200_min_bias_idnni_nx_201.root

#IDNNI PARAMS
xmax 12
ymax 12
nx   241
ny   241

gaussian_smearing_sigma 0.1

upper_mult_proxy_cut  1000.00
lower_mult_proxy_cut  257.995




#rapidity extension
etamax 10
neta   101
matter_eta_plateau 1.3
matter_eta_fall  1.5
matter_tilt_param  2.5

baryon_rapidity_profile_eta_peak  3.5
baryon_rapidity_profile_sigma_eta_plus  0.1
baryon_rapidity_profile_sigma_eta_minus  2.0



# :: END :: #
