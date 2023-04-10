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

# To run this example, we need to set the working directory to the folder
# where the https://github.com/imarinfr/klo repository has been checked out.
#setwd(<CHECKOUT DIRECTORY>/example)

# The first step is to load and reformat each data set into three columns
# of LMS values lms1 and lms2
library(klo)
library(R.matlab)
lms1 <- readMat("../data/lms_sete_fontes_1320.mat")$lms1
lms2 <- readMat("../data/lms_sete_fontes_1321.mat")$lms2
nr <- dim(lms1)[1]
nc <- dim(lms1)[2]
nw <- dim(lms1)[3]
dim(lms1) <- c(nr * nc, nw)
dim(lms2) <- c(nr * nc, nw)

# Here are the first few rows of lms1 and lms2. Values are less than unity
# because of the way LMS values are normalized.
head(lms1)
head(lms2)

# The initial Gaussian estimate of differential entropy for the first image
# is obtained thus:
entg(lms1)

# These values are negative because the data values are less than unity.
# The KL estimates of differential entropy without and then with offset are
# obtained similarly:
entkl(lms1, "kl")$ent
entkl(lms1, "klo")$ent

# The initial Gaussian estimate of mutual information between the two
# images is obtained thus:
mig(lms1, lms2)

# The estimate is non-negative, as expected. The KL estimates of mutual
# information without and then with offset are obtained similarly:
mikl(lms1, lms2, "kl")$mi
mikl(lms1, lms2, "klo")$mi

# Notice the progressive increase in the size of the mutual information
# estimates with the KL estimator and the KLo estimator

###########################################################################
# (c) Ivan Marin-Franch and David H Foster 12-May-21