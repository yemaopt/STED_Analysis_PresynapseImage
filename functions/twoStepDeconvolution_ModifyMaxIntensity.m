function [deconvolvedImData, enhancedPSF] = twoStepDeconvolution_ModifyMaxIntensity(imData,PSF,numberIterations)
%% twoStepDeconvolution
% imagePath is the path to the folder containing the image file.
% imageName is the name of the image to deconvolve
% psfPath is the path to the PSF file - e.g., /project/cellbiology/Dean_lab/shared/psfs/ctASLM2-510nm.tif
% numberIterations is usually set to 10.
%
% Written by Bo-Jui Chang, 2019.  Verified on Matlab/2019a.

% Threshold PSF by bottom 5% & Normalize
intensityDistribution = sort(PSF(:));
PSFbackground = mean(intensityDistribution(1:size(PSF(:))/20));
disp(['The Background Intensity for the PSF is ' num2str(PSFbackground)]);
PSF=abs(PSF-PSFbackground);

[ny,nx]=size(imData);

% Deconvolve the PSF to get a better estimate of the real PSF.
% Load the data.

paddedImData=padarray(single(imData),[20 20],'symmetric'); 
disp('Deconvolving Data');
[~,enhancedPSF]=deconvblind(paddedImData,PSF,numberIterations);

%% Deconvolve the Data With The Improved PSF Estimate
disp('Deconvolving Data with Enhanced PSF');
[deconvolvedImData,~]=deconvblind(paddedImData,enhancedPSF,numberIterations);
deconvolvedImData=deconvolvedImData(21:20+ny,21:20+nx);
deconvolvedImData=deconvolvedImData./max(deconvolvedImData(:));
deconvolvedImData=uint16(deconvolvedImData*2^16);


