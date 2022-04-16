function [] = write2nifti( vol, data, outfile, dims )
% 
% [] = write2nifti( vol, data, outfile[, dims] )
% 

if nargin<4
    dims = size(data);
end

vol.img = data;
vol.hdr.dime.dim(1) = length(dims);
vol.hdr.dime.dim(2:2+length(dims)-1) = dims;
dmax = prctile(data(data>0),95);
if length(dims)<4
    vol.hdr.dime.cal_min = 0;
    vol.hdr.dime.cal_max = dmax;
    vol.hdr.dime.glmin = 0;
    vol.hdr.dime.glmax = dmax;
end
save_untouch_nii(vol,outfile);

end