# Script to obtain sample estimates of differential entropy and mutual
# information from trivariate data. Example 1 is the same as in the article
# apart from image subsampling. Example 2 focuses on processing scene data
# to illustrate effects of different illuminations on information.

# set the working directory
setwd("~/06.optocom/03.software/klo_v1/R/")

# load the necessary libraries
library(FNN)      # for nearest-neighbor search
library(expm)     # for the matrix square root
library(R.matlab) # to load MatLab files in R

# source the codes with need to compute the estimates
source("mikl.r")
source("rgbcurves.r")

# EXAMPLE 1
# Load 1005 x 1306 x 3 scene data coded as LMS (long-, medium-, and
# short-wavelength) cone photoreceptor responses of human eye at each pixel
lms1 <- readMat("../data/lms_sete_fontes_1333.mat")$lms1
lms2 <- readMat("../data/lms_sete_fontes_1335.mat")$lms2

# Reshape LMS scene data into 1312530 x 3 arrays
nr <- dim(lms1)[1]
nc <- dim(lms1)[2]
nw <- dim(lms1)[3]
dim(lms1) <- c(nr * nc, nw)
dim(lms2) <- c(nr * nc, nw)

# Get Gaussian and KL estimates of differential entropy without and
# with offset correction
print(entg(cbind(lms1, lms2)))
print(entkl(cbind(lms1, lms2), "kl"))
print(entkl(cbind(lms1, lms2), "klo"))

# Get Gaussian and KL estimates of mutual information without and
# with the offset correction
print(mig(lms1, lms2))
print(mikl(lms1, lms2, "kl"))
print(mikl(lms1, lms2, "klo"))

###########################################################################
# EXAMPLE 2
# Now try similar calculation but using hyperspectral reflectance images of
# scene under different illuminants and viewed by camera with different
# sensor spectra or by human eye. For more details see following:
# https://personalpages.manchester.ac.uk/staff/david.foster/Tutorial_Colour_Information/Tutorial_Color_Information.html

# Load 255 x 335 x 33 scene spectral reflectance data (33 wavelengths)
reflectances <- readMat("../data/ref4_scene5.mat")$reflectances
nr <- dim(reflectances)[1]
nc <- dim(reflectances)[2]
nw <- dim(reflectances)[3]

# Load daylight illuminants of correlated color temperature 25000 K and
# 4000 K; other daylight illuminant combinations can be tested for their
# effect on mutual information, e.g. 6500 K
illum_25000 <- readMat("../data/illum_25000.mat")$illum.25000
illum_4000 <-readMat("../data/illum_4000.mat")$illum.4000

# Convert to 255 x 335 x 33 scene spectral radiances for each illuminant
# (wavelength for-loop is used for clarity, not efficiency)
radiances_25000 <- array(rep(NA, nr * nc * nw), dim = c(nr, nc, nw))
radiances_4000  <- array(rep(NA, nr * nc * nw), dim = c(nr, nc, nw))
for(i in 1:33) {
  radiances_25000[,,i] <- reflectances[,,i] * illum_25000[i]
  radiances_4000[,,i]  <- reflectances[,,i] * illum_4000[i]
}

# Reshape scene data into 85425 x 3 arrays
dim(radiances_25000) <- c(nr * nc, nw)
dim(radiances_4000)  <- c(nr * nc, nw)

# Get 33 x 4 RGB sensor curves for chosen sensor, here 'agilent'
# Remove first column wavelength vector as implicit in data structure
rgbsens <- rgbcurves("agilent")[,-1]
rgb_25000 <- radiances_25000 %*% rgbsens
rgb_4000  <- radiances_4000 %*% rgbsens

# Gaussian and KL estimates of differential entropy without and
# with the offset correction
print(entg(cbind(rgb_25000, rgb_4000)))
print(entkl(cbind(rgb_25000, rgb_4000), 'kl'))
print(entkl(cbind(rgb_25000, rgb_4000), 'klo'))

# Gaussian and KL estimates of mutual information without and
# with the offset correction
print(mig(rgb_25000, rgb_4000))
print(mikl(rgb_25000, rgb_4000, 'kl'))
print(mikl(rgb_25000, rgb_4000, 'klo'))

###########################################################################
# (c) Ivan Marin-Franch and David H Foster 10-Jul-20