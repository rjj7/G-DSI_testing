
%% Setup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isdeployed
    disp('--RUNNING DEPLOYED APP--');
else
    addpath('/autofs/space/hemera_001/users/rjjones/code/toolboxes/NIfTI_20140122');
    addpath(genpath('/autofs/space/nyx_002/users/rjones/Bay8/code/gdsi/version_020122'));
end
% clear; close all; clc;
verbose = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Hard-coded variables ---
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%diffusivity of free water in tissue at 37C
%   (this is to set MDD_water, for F_mtx)
D_water = 2.5e-3;

%use spherical vertices from DSI studio (8-fold tesselation)
%   (also used to create F_mtx)
Use_DSI_Studio_dirs     = true;

%if 2+ lowb vols, concat and place in front of 4D dwi stack.
%   (and do the same for bvals/bvecs - 
%    this is needed to create qvecs, then precompute F_mtx)
concatlowb              = true;     %concat + avg lowb volumes


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- input data ---
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   inputs/paths
dtroot      = '/autofs/space/nyx_002/users/rjones/Bay8/CS_DSI_10_01_2021';

fdwi        = [dtroot filesep 'preproc' filesep 'dwi_eddy.nii.gz'];
fmask       = [dtroot filesep 'preproc' filesep 'brain_mask.nii.gz'];
fbvec       = [dtroot filesep 'preproc' filesep 'dwi_eddy.eddy_rotated_bvecs'];
fbval       = [dtroot filesep 'preproc' filesep 'bvals'];
flowb       = [dtroot filesep 'preproc' filesep 'lowb_eddy_mean.nii.gz'];
outdir      = [dtroot filesep 'preproc' filesep 'analysis' filesep 'gdsi' filesep 'results-testing'];
procdir     = [dtroot filesep 'preproc' filesep 'preproc'];
acqspecs    = [dtroot filesep 'preproc' filesep 'analysis' filesep 'gdsi' filesep 'acqspecs.txt'];

big_delta   = 25.5 * 1e-3;    %ms
small_delta = 10.8 * 1e-3;    %ms




%% Load data               

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Load diffusion data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if verbose
    disp('Loading subject data...');
end
% [ DATAIN ] = load_subject_data( INFO, verbose ); %opts.loaddwi=true
DATAIN = [];

%%%% DWI (4D VOL)
if verbose
    disp('Loading dwi...');
end
vol = load_untouch_nii(fdwi);
DATAIN.dwi = double(vol.img);
% DATAIN.lowb0 = DATAIN.dwi(:,:,:,1);

%%%% BRAIN MASK (3D VOL)
if verbose
    disp('Loading mask...');
end
vol = load_untouch_nii(fmask);
DATAIN.mask = double(vol.img);

%%%% BVEC/BVAL
if verbose
    disp('Loading bvecs/bvals...');
end
DATAIN.bvecs = dlmread(fbvec);
DATAIN.bvals = dlmread(fbval);

% Check that dims of bvecs/bvals/dwi match up
validate_diff_data( DATAIN );
if verbose
    disp(' -- Diffusion data validated --');
end

%%%% Find b=0 indices
index_b0 = DATAIN.bvals==0;
DATAIN.index_b0 = index_b0;

%%%% MEAN LOWB VOLUME
% DATAIN.lowb = mean(DATAIN.dwi(:,:,:,index_b0),4);

%%% Remove duplicate lowb volumes
if verbose
    disp('Creating qvec...');
end
if size(DATAIN.bvals,2)>size(DATAIN.bvals,1)
    DATAIN.bvals = DATAIN.bvals';
end
if size(DATAIN.bvecs,2)>size(DATAIN.bvecs,1)
    DATAIN.bvecs = DATAIN.bvecs';
end
if concatlowb
    DATAIN.bvals_orig = DATAIN.bvals;
    DATAIN.bvecs_orig = DATAIN.bvecs;
    DATAIN.bvals = DATAIN.bvals(~index_b0);
    DATAIN.bvals = cat(1,0,DATAIN.bvals);
    DATAIN.bvecs = DATAIN.bvecs(~index_b0,:);
    DATAIN.bvecs = cat(1,[0 0 0],DATAIN.bvecs);
end

%%% Create qvec
qvec = repmat(sqrt( 6 * D_water * DATAIN.bvals ), [1, 3]) .* DATAIN.bvecs;
nqvec = round(qvec * (5/round(max(qvec(:)))));
mag_qvec = vecnorm(nqvec,2,2);

DATAIN.qvec = qvec;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Prepare for GDSI recon
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% Calculate mean diffusion distance of water (micron)
mdd_water = sqrt(6 * D_water * (big_delta-(small_delta/3)) ) * 1000;
% mdd_water = 10;
fprintf('   Diffusivity \t = %g mm^2/s\n', D_water);
fprintf('   MDD water\t = %.2f micron\n', mdd_water);


%% Generate vertices on sphere
%%% Define EAP directions on sphere
if verbose
    disp('Creating pdf_dirs...');
end
if Use_DSI_Studio_dirs
    % Use DSI Studio dirs
    load('odf8.mat');
    diff_disp_sphere.vertices = odf_vertices';
    diff_disp_sphere.faces = odf_faces';
    nb_pdf_dirs_full           = length(diff_disp_sphere.vertices);
    nb_pdf_dirs                = nb_pdf_dirs_full/2;
    pdf_dirs              = diff_disp_sphere.vertices(1:nb_pdf_dirs,:);
    
    odf_verticesI = odf_vertices';
    odf_facesI = odf_faces';
else
    % Use matlab sphere function 
%     [ diff_disp_sphere ]  = create_sphere( 25 );
    [ diff_disp_sphere ]  = create_sphere_nneg( 25 );
    nb_pdf_dirs_full      = length(diff_disp_sphere.vertices);
    nb_pdf_dirs           = nb_pdf_dirs_full;
    pdf_dirs              = diff_disp_sphere.vertices(1:nb_pdf_dirs,:);

end


%% Generate radial displacement vector
%%% Define displacement vector
if verbose
    disp('Creating displacement vec...');
end
% Define radial points to evaluate EAP at
nr = 100; %number of points
rs = linspace(0, 1, nr)';  
dispvec  = linspace(0, mdd_water, length(rs));
r_step_size           = dispvec(2)-dispvec(1);
dispvec2 = dispvec.^2;
dispvec4 = dispvec.^4;
        

%% Precompute Fmtx
% Original Fmtx (RxNxD)
tic1 = tic;
[ Fmtx ] = precompute_Fmtx( qvec, pdf_dirs, rs );          
toc1 = toc(tic1);
fprintf(' time to compute Fmtx = %g sec\n',toc1);

%New concat Fmtx (RDxN)
tic2 = tic;
N = length(qvec);
D = size(pdf_dirs,1);
Fmtx2 = zeros(nr*D,N);
for pdfind=1:nr
    pstart = (pdfind-1)*D + 1;
    pend = pstart + D - 1;
    
    Fmtx2(pstart:pend,:) = reshape(Fmtx(pdfind,:,:),D,N );
end
toc2 = toc(tic2);
fprintf(' time to compute Fmtx2 = %g sec\n',toc2);

disp(' ');

% xrange = 52:70;
% yrange = 38:56; 
% zrange = 35;

xrange = 1:140;
yrange = 1:140;
zrange = 35;


% [ ~, res1 ] = debug_GDSI_recon_v1( DATAIN, xrange, yrange, zrange, Fmtx );
[ ~, res2 ] = debug_GDSI_recon_v2( DATAIN, xrange, yrange, zrange, Fmtx2 );


disp(' ');



