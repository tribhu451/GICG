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
SNN 2760.0


#collision species & impact parameter
#[Info] species can be (Au,U,Pb,p) 
projectile Pb
target Pb
bmin  7.0001
bmax  7.0002

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


EMsigma 0.023


# :: END :: #






