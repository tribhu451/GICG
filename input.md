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
projectile Au_trento
target Au_trento
bmin  7.0001
bmax  7.0001


gamma_fluctuation_k 1.4

#Grid, smearing and event averaged profile parameters
xmax 12
ymax 12
nx   241
ny   241
gaussian_smearing_sigma 0.4

write_profile 1


upper_mult_proxy_cut  10000000.00
lower_mult_proxy_cut  00000000.00


# :: END :: #






