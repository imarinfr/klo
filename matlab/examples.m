% Script to obtain sample estimates of differential entropy and mutual
% information from trivariate data. Example 1 is the same as in the article
% apart from image subsampling. Example 2 focuses on processing scene data
% to illustrate effects of different illuminations on information.

% EXAMPLE 1
% Load 1005 x 1306 x 3 scene data coded as LMS (long-, medium-, and
% short-wavelength) cone photoreceptor responses of human eye at each pixel
load ../data/lms_sete_fontes_1333.mat; % returns lms1 at 13:33 h
load ../data/lms_sete_fontes_1335.mat; % returns lms2 at 13:35 h

% Reshape LMS scene data into 1312530 x 3 arrays
[nr, nc, nw] = size(lms1);
lms1 = reshape(lms1, nr * nc, nw);
lms2 = reshape(lms2, nr * nc, nw);

% Get Gaussian and KL estimates of differential entropy without and
% with offset correction
entg([lms1 lms2])
entkl([lms1 lms2], 'kl')
entkl([lms1 lms2], 'klo')

% Get Gaussian and KL estimates of mutual information without and
% with the offset correction
mig(lms1, lms2)
mikl(lms1, lms2, 'kl')
mikl(lms1, lms2, 'klo')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EXAMPLE 2
% Now try similar calculation but using hyperspectral reflectance images of
% scene under different illuminants and viewed by camera with different
% sensor spectra or by human eye. For more details see following:
% https://personalpages.manchester.ac.uk/staff/david.foster/Tutorial_Colour_Information/Tutorial_Color_Information.html

% Load 255 x 335 x 33 scene spectral reflectance data (33 wavelengths)
load ../data/ref4_scene5.mat;
[nr, nc, nw] = size(reflectances); % returns reflectances

% Load daylight illuminants of correlated color temperature 25000 K and
% 4000 K; other daylight illuminant combinations can be tested for their
% effect on mutual information, e.g. 6500 K
load ../data/illum_25000.mat;   % returns illum_2500
load ../data/illum_4000.mat;    % returns illum_4000

% Convert to 255 x 335 x 33 scene spectral radiances for each illuminant
% (wavelength for-loop is used for clarity, not efficiency)
radiances_25000 = zeros(nr, nc, nw);
radiances_4000 = zeros(nr, nc, nw);
for i = 1:33
  radiances_25000(:,:,i) = reflectances(:,:,i) * illum_25000(i);
  radiances_4000(:,:,i) = reflectances(:,:,i) * illum_4000(i);
end

% Reshape scene data into 85425 x 3 arrays
radiances_25000 = reshape(radiances_25000, nr * nc, nw);
radiances_4000 = reshape(radiances_4000, nr * nc, nw);

% Get 33 x 4 RGB sensor curves for chosen sensor, here 'agilent'
rgbsens = rgbcurves('agilent'); 
% Remove first column wavelength vector as implicit in data structure
rgbsens(:,1) = []; 

% Get 85425 x 3 RGB values for reshaped radiance image 
rgb_25000 = radiances_25000 * rgbsens;   
rgb_4000 = radiances_4000 * rgbsens;

% Get Gaussian and KL estimates of differential entropy without and
% with the offset correction
entg([rgb_25000 rgb_4000])
entkl([rgb_25000 rgb_4000], 'kl')
entkl([rgb_25000 rgb_4000], 'klo')

% Get Gaussian and KL estimates of mutual information without and
% with the offset correction
mig(rgb_25000, rgb_4000)
mikl(rgb_25000, rgb_4000, 'kl')
mikl(rgb_25000, rgb_4000, 'klo')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Ivan Marin-Franch and David H Foster 10-Jul-20