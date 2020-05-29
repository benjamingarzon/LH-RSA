% Scripts to remove residual B1 bias from T1 maps calculated with the
% MP2RAGE sequence
% Correction of T1 estimations as sugested on:
% Marques, J.P., Gruetter, R., 2013. New Developments and Applications of the MP2RAGE Sequence - Focusing the Contrast and High Spatial Resolution R1 Mapping. PLoS ONE 8. doi:10.1371/journal.pone.0069294
%
% in this script it is assumed that the B1 maps have been coregistered to the
% space of the MP2RAGE image and that they now have the B1 has in the
% process been interpolated to the same resolution.
% 

function estimateT1(MP2name, T1outname, varargin);
addpath('/home/benjamin.garzon/Software/LeftHand/process/mprageconvert/Nifti_tools')
addpath(genpath('.'))
if nargin==1
    B1name = varargin{1};
else
    B1name = "";
end
display(T1outname)
 %MP2name = '/home/benjamin.garzon/Data/LeftHand/Lund1/vbm/sub-lue1101/ses-1/T1wtotemplate_brain.nii.gz';
 %T1outname = '/home/benjamin.garzon/Data/LeftHand/Lund1/vbm/sub-lue1101/ses-1/T1map.nii.gz';
 invEFF = 0.8;%0.96;

    %% MP2RAGE protocol info and loading the MP2RAGE dataset 
    
    MP2RAGE.B0=7;           % in Tesla
    MP2RAGE.TR=5;           % MP2RAGE TR in seconds 
    MP2RAGE.TRFLASH=[6.8e-3]; % TR of the GRE readout
    MP2RAGE.TIs=[0.911 2.761]; %inversion times - time betweenC middle of refocusing pulse and excitatoin of the k-space center encoding
    MP2RAGE.FlipDegrees=[5 3];% Flip angle of the two readouts in degrees
    MP2RAGE.filename = MP2name; % file
    
    % load the MP2RAGE data
    MP2RAGEimg=load_untouch_nii(MP2RAGE.filename);
    MP2RAGE.NZslices=[ size(MP2RAGEimg.img, 3) ];% Slices Per Slab * [PartialFourierInSlice-0.5  0.5]

    % check the properties of this MP2RAGE protocol... this happens to be a
    % very B1 insensitive protocol
    plotMP2RAGEproperties(MP2RAGE)
    
%% if another technique was used to obtain the relative B1 maps 
%  (1 means correctly calibrated B1)
    if B1name ~= "" 
      B1=load_untouch_nii(B1name);      
      B1.img = double(B1.img)/100;
    else
      B1 = MP2RAGEimg;
      B1.img = double(1.*(B1.img>0));
    end
    brain = MP2RAGEimg;
    MP2RAGEimg.img =  MP2RAGEimg.img - 0.5; % correct the offset
    [ T1corrected, MP2RAGEcorr] = T1B1correctpackageTFL(B1, MP2RAGEimg,[], MP2RAGE, brain, invEFF);
    %T1corrected.img(~brain.img) = 0;
    % saving the data
%    save_untouch_nii(MP2RAGEcorr,MP2outname) 
    save_untouch_nii(T1corrected, T1outname)

end