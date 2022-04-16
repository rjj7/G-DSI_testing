function [fslut] = load_fs_lut(lutfile)
%
% [fslut] = load_fs_lut(lutfile)
% 
% Loads FreeSurfer LUT for aseg+aparc.
% 
% If no arguments given, loads default file.
% 
% Outputs an array fslut with 4 columns.
%  1st column is the label numbers
%  2nd-4th columns are the R, G and B channel intensities 
% 

if nargin<1
    lutfile = '/autofs/space/hemera_002/users/rjjones/gdap/cfg/fs/matlab_FS_LUT_table.txt';
end

temp = readtable(lutfile); %load modified lut file
fslut.labelnb = table2array(temp(:,1));
fslut.labelname = table2cell(temp(:,2));
fslut.labelrgb = table2array(temp(:,3:5))/255; 
% fslut = table2array(fslut(:,[1,3:5])); %grab label and rgb cols
% fslut(:,2:4) = fslut(:,2:4)/255; %to normalize from max=255 (uint8) to max=1

end