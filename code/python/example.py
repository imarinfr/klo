# This example shows how to compute the Gaussian approximation and the KL
# and KLo estimates of differential entropy and of mutual information
# between two sets of trivariate data. Each set represents a color image of
# a scene recorded at successive instants, about 1 minute apart, shown in
# Figure 1 of the manuscript (other images are available with the
# software). The image values were expressed not as conventional RGB
# triplets but as LMS triplets, corresponding to activities in the long-,
# medium-, and short-wavelength-sensitive cone photoreceptors of the eye.
# Each data set is stored as a 507x657x3 array, where the first two
# dimensions index pixel coordinates and the third dimension indexes LMS
# values.

# Import functions necessary to run entkl and mikl
from numpy import empty, reshape, concatenate
from scipy.io import loadmat
from mikl import entg, mig, entkl, mikl

# To run this example, we need to set the working directory to the
# subfolder code/python/ in the folder where the
# https://github.com/imarinfr/klo repository has been checked out.

# The first step is to load and reformat each data set into three columns
# of LMS values lms1 and lms2
lms1 = loadmat("../../data/lms_sete_fontes_1320.mat")["lms1"]
lms2 = loadmat("../../data/lms_sete_fontes_1321.mat")["lms2"]
nr, nc, nw = lms1.shape
lms1 = reshape(lms1, (nr * nc, nw), order="F")
lms2 = reshape(lms2,  (nr * nc, nw), order="F")

# Here are the first few rows of lms1 and lms2. Values are less than unity
# because of the way LMS values are normalized.
print(lms1[:6,:])
print(lms2[:6,:])

# The initial Gaussian estimate of differential entropy for the first image
# is obtained thus:
print(entg(lms1))

# These values are negative because the data values are less than unity.
# The KL estimates of differential entropy without and then with offset are
# obtained similarly:
print(entkl(lms1, "kl"))
print(entkl(lms1, "klo"))

# The initial Gaussian estimate of mutual information between the two
# images is obtained thus:
print(mig(lms1, lms2))

# The estimate is non-negative, as expected. The KL estimates of mutual
# information without and then with offset are obtained similarly:
print(mikl(lms1, lms2, "kl"))
print(mikl(lms1, lms2, "klo"))

# Notice the progressive increase in the size of the mutual information
# estimates with the KL estimator and the KLo estimator

###########################################################################
# (c) Ivan Marin-Franch and David H Foster 12-May-21