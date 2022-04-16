function validate_diff_data( DATA )
% 
% validate_diff_data( DATA )
% 
% Checks that # DWIs matches lines of bvec, bval;
% error if mismatch
%

nframes = size(DATA.dwi,4);
nbvals = length(DATA.bvals);
nbvecs = length(DATA.bvecs);

if nbvecs ~= nbvals 
    error('diffMismatch:BveBva','\nDims mismatch: bvecs + bvals\n');
end
if nbvecs ~= nframes 
    error('diffMismatch:BveDwi','\nDims mismatch: bvecs + dwis\n');
end
if nbvals ~= nframes
    error('diffMismatch:BvaDwi','\nDims mismatch: bvals + dwis\n');
end

end