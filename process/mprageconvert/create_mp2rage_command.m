% Create MP2RAGE image volume as .nii.gz from one real and imaginary
% MP2RAGE data sets. Zipped 4D NIfTI files are needed. For scaling from
% stored values to floating point values one DICOM from the acquisition
% need to be specified.
% MP2RAGE: Marques et al. NeuroImage 49 (2010) 1271-1281 MP2RAGE, a self
% bias-field corrected sequence for improved segmentation and T1-mapping at
% high field.


function create_mp2rage_command(mypath, re_file, im_file, jsonfile);

re_path = mypath;
im_path = mypath;

% Load .nii.gz data
re_nii = load_untouch_nii(strcat(re_path,re_file));
im_nii = load_untouch_nii(strcat(im_path,im_file));

% Get image data
re_ti1_sv = double(re_nii.img(:,:,:,1));
re_ti2_sv = double(re_nii.img(:,:,:,2));
im_ti1_sv = double(im_nii.img(:,:,:,1));
im_ti2_sv = double(im_nii.img(:,:,:,2));

% Get scale intercept and scale slope from DICOM header for scaling from
% stored values to floating point values
val = jsondecode(fileread(strcat(mypath,jsonfile)));
ri = val.PhilipsRescaleIntercept; %;dcm_header.Private_2005_100d;
rs = val.PhilipsRescaleSlope; %dcm_header.Private_2005_100e;
ss = val.PhilipsScaleSlope; 
si = -ri/rs;

%Scale from sv to fp
re_ti1 = double((re_ti1_sv-si)/ss);
re_ti2 = double((re_ti2_sv-si)/ss);
im_ti1 = double((im_ti1_sv-si)/ss);
im_ti2 = double((im_ti2_sv-si)/ss);  

% Create complex image matrices
ti1 = re_ti1+1i.*im_ti1;
ti2 = re_ti2+1i.*im_ti2;

% Marquez 2010 Eq.(3)
MP2RAGE = real(conj(ti1).*ti2./(abs(ti1).^2+abs(ti2).^2));

% Save as NIfTI file, change data type, and remove NIfTI scaling
re_nii.img = MP2RAGE;
re_nii.hdr.dime.dim = [3 re_nii.hdr.dime.dim(2) re_nii.hdr.dime.dim(3) re_nii.hdr.dime.dim(4) 1 1 1 1];
re_nii.hdr.dime.scl_slope = 1;
re_nii.hdr.dime.scl_inter = 0;
re_nii.hdr.dime.datatype = 16;    %real
re_nii.hdr.dime.bitpix = 32;      %32bit
save_untouch_nii(re_nii, strcat(re_path,'MP2RAGE.nii.gz'))


