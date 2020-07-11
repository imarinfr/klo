# Script to obtain sample estimates of differential entropy and mutual
# information from trivariate data. Example 1 is the same as in the article
# apart from image subsampling. Example 2 focuses on processing scene data
# to illustrate effects of different illuminations on information.

# Import functions necessary to run entkl and mikl
from numpy import empty, reshape, concatenate
from scipy.io import loadmat
from rgbcurves import rgbcurves
from mikl import entg, mig, entkl, mikl

# EXAMPLE 1
# Load 1005 x 1306 x 3 scene data coded as LMS (long-, medium-, and
# short-wavelength) cone photoreceptor responses of human eye at each pixel
lms1 = loadmat("../data/lms_sete_fontes_1333.mat")["lms1"]
lms2 = loadmat("../data/lms_sete_fontes_1335.mat")["lms2"]
nr, nc, nw = lms1.shape

# Reshape LMS scene data into 1312530 x 3 arrays
lms1 = reshape(lms1, (nr * nc, nw), order="F")
lms2 = reshape(lms2,  (nr * nc, nw), order="F")

# Get Gaussian and KL estimates of differential entropy without and
# with offset correction
print(entg(concatenate((lms1, lms2), axis=1)))
print(entkl(concatenate((lms1, lms2), axis=1), "kl"))
print(entkl(concatenate((lms1, lms2), axis=1), "klo"))

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
reflectances = loadmat("../data/ref4_scene5.mat")["reflectances"]
nr, nc, nw = reflectances.shape

# Load daylight illuminants of correlated color temperature 25000 K and
# 4000 K; other daylight illuminant combinations can be tested for their
# effect on mutual information, e.g. 6500 K
illum_25000 = loadmat("../data/illum_25000.mat")["illum_25000"]
illum_4000 = loadmat("../data/illum_4000.mat")["illum_4000"]

# Convert to 255 x 335 x 33 scene spectral radiances for each illuminant
# (wavelength for-loop is used for clarity, not efficiency)
radiances_25000 = empty((nr, nc, nw))
radiances_4000 = empty((nr, nc, nw))
for i in range(33):
    radiances_25000[:, :, i] = reflectances[:, :, i] * illum_25000[i]
    radiances_4000[:, :, i] = reflectances[:, :, i] * illum_4000[i]

# Reshape scene data into 85425 x 3 arrays
radiances_25000 = reshape(radiances_25000, (nr * nc, nw), order="F")
radiances_4000 = reshape(radiances_4000,  (nr * nc, nw), order="F")


# Get 33 x 4 RGB sensor curves for chosen sensor, here 'agilent'
# Remove first column wavelength vector as implicit in data structure
rgbsens = rgbcurves("agilent")[:,1:]
rgb_25000 = radiances_25000.dot(rgbsens)
rgb_4000 = radiances_4000.dot(rgbsens)

# Gaussian and KL estimates of differential entropy without and
# with the offset correction
print(entg(concatenate((rgb_25000, rgb_4000), axis=1)))
print(entkl(concatenate((rgb_25000, rgb_4000), axis=1), 'kl'))
print(entkl(concatenate((rgb_25000, rgb_4000), axis=1), 'klo'))

# Gaussian and KL estimates of mutual information without and
# with the offset correction
print(mig(rgb_25000, rgb_4000))
print(mikl(rgb_25000, rgb_4000, 'kl'))
print(mikl(rgb_25000, rgb_4000, 'klo'))

###########################################################################
# (c) Ivan Marin-Franch and David H Foster 10-Jul-20