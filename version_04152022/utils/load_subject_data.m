function [ DATA ] = load_subject_data( INFO, verbose, opts )
% [ DATA ] = load_subject_data( INFO, verbose, opts )

if ~isdeployed
    addpath('/autofs/space/hemera_001/users/rjjones/code/toolboxes/NIfTI_20140122');
end

if nargin<2
    verbose = 0;
end

if nargin<3
    opts = [];
    opts.loaddwi=true;
end


%%%%
DATA = [];


%%%% DWI
if isfield(opts,'loaddwi') && opts.loaddwi
    if verbose
        disp('Loading dwi...');
    end
    vol = load_untouch_nii(INFO.fdwi);
    DATA.dwi = double(vol.img);
    DATA.lowb = DATA.dwi(:,:,:,1);
else
    DATA.dwi = [];
    if verbose
        disp(' - Not loading dwi');
    end
    if isfield(INFO,'flowb') && exist(INFO.flowb,'file')
        disp('  (loading lowb instead!) ');
        vol = load_untouch_nii(INFO.flowb);
        DATA.lowb = double(vol.img);
    end
end


%%%% BRAIN MASK
if verbose
    disp('Loading mask...');
end
vol = load_untouch_nii(INFO.fmask);
DATA.mask = double(vol.img);


%%%% BVEC/BVAL
if verbose
    disp('Loading bvecs/bvals...');
end
DATA.bvecs = dlmread(INFO.fbvec);
DATA.bvals = dlmread(INFO.fbval);


end



