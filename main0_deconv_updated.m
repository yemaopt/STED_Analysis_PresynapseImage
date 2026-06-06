% Written by Ye Ma, Department of Biomedical Engineering, Johns Hopkins University, 2022

close all
clear
clc
addpath ./functions/bfmatlab
addpath ./functions

imageDirectory = uigetdir('', 'Please select the folder containing raw prp images');
imglist=dir(fullfile(imageDirectory,sprintf('*.prp')));

numberIterations=10;
outputDirectory=fullfile(imageDirectory,strcat('decon_',num2str(numberIterations))); mkdir(outputDirectory);

[filename_635, pathname_635, ~] = uigetfile('*.tif', 'Please select the PSF for 635 channel');
[filename_594, pathname_594, ~] = uigetfile('*.tif', 'Please select the PSF for 594 channel');

PSF_635 = double(tiffRead(fullfile(pathname_635,filename_635)));
PSF_594 = double(tiffRead(fullfile(pathname_594,filename_594)));

for ii = 1:length(imglist)
    imageName = imglist(ii).name;
    data = bfopen([imglist(ii).folder,'/',imglist(ii).name]);
    img_ch1=cell2mat(data{1,1}(1));
    img_ch2=cell2mat(data{2,1}(1));
    img_635=img_ch1;img_635=imgaussfilt(img_635,1.2);
    img_594=img_ch2;img_594=imgaussfilt(img_594,1.2);
%     [img_635_deconv,enhancedPSF_635]=twoStepDeconvolution_ModifyMaxIntensity(img_635,PSF_635,numberIterations);
%     [img_594_deconv,enhancedPSF_594]=twoStepDeconvolution_ModifyMaxIntensity(img_594,PSF_594,numberIterations);
    
    [img_635_deconv,enhancedPSF_635]=twoStepDeconvolution(img_635,PSF_635,numberIterations);
    [img_594_deconv,enhancedPSF_594]=twoStepDeconvolution(img_594,PSF_594,numberIterations);

    img_composite_deconv(:,:,1)=img_635_deconv;
    img_composite_deconv(:,:,2)=img_594_deconv;
    img_composite_deconv=uint16(img_composite_deconv);
    write3Dtiff(img_composite_deconv, fullfile(outputDirectory,strcat(imageName,'_DualBlindDecon_',num2str(numberIterations),'.tif')))
    tiffWrite(uint16(enhancedPSF_635./max(enhancedPSF_635(:)).*2^16),fullfile(outputDirectory,'enhancedPSF_635.tif'));
    tiffWrite(uint16(enhancedPSF_594./max(enhancedPSF_594(:)).*2^16),fullfile(outputDirectory,'enhancedPSF_594.tif'));
end
