function [ DATAOUT ] = script_GDSI_recon_HCP_3sh( verbose ) %( fconfig, verbose )
% 
% GDSI reconstruction of RSI data
%
% original script by (c) Qiyuan Tian, Stanford RSL
% Modified by RJ, Feb 2021 - Feb 2022
%

% % --- EXAMPLE ---
% % fconfig = '/autofs/space/nyx_002/users/rjones/Bay8/code/gdsi/version_020122/hcp_mgh.config';
% fconfig = '/autofs/space/nyx_002/users/rjones/Bay8/CS_DSI_10_01_2021/analysis/gdsi/dsi.config';
% verbose = true;


%% Setup

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isdeployed
    disp('--RUNNING DEPLOYED APP--');
else
    addpath('/autofs/space/hemera_001/users/rjjones/code/toolboxes/NIfTI_20140122');
    addpath(genpath('/autofs/space/nyx_002/users/rjones/Bay8/code/gdsi/version_020122'));
end


% clear; close all; clc;

if nargin<1, verbose = true; end



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


%other parameters
Remove_HighBValue       = false;    %discard high-bvals
ClipRmatrix             = false;    %clip negative values of propagator
Norm2PDF                = false;    %normalize EAP so sum(EAP(:))=1
ReconODFs               = false;    %recon dODFs
ReconZDPs               = false;    %recon zero-displacement maps (RTOP, RTAP, RTPP)
ReconMDs                = false;    %recon mean displacements (1st, 2nd, 4th moments)
MakePlots.odf           = false;    %make plot of ODF
MakePlots.rgbprofiles   = false;    %make direction-encoded plot of all 1D EAP profiles
StoreRmatrix            = true;     %store full EAP in each voxel (ndirs-by-nr mtx)


%   inputs/paths
dtroot      = '/autofs/space/nyx_002/users/rjones/Bay8/CS_DSI_10_01_2021';

fdwi        = [dtroot filesep 'preproc' filesep 'nufft/3sh' filesep 'data.nii.gz'];
fmask       = [dtroot filesep 'preproc' filesep 'nufft/3sh' filesep 'brain_mask.nii.gz'];
fbvec       = [dtroot filesep 'preproc' filesep 'nufft/3sh' filesep 'bvecs'];
fbval       = [dtroot filesep 'preproc' filesep 'nufft/3sh' filesep 'bvals'];
% flowb       = [dtroot filesep 'preproc' filesep 'lowb_eddy_mean.nii.gz'];
outdir      = [dtroot filesep 'preproc' filesep 'analysis' filesep 'gdsi' filesep 'results-testing'];
% procdir     = [dtroot filesep 'preproc' filesep 'preproc'];
% acqspecs    = [dtroot filesep 'preproc' filesep 'analysis' filesep 'gdsi' filesep 'acqspecs.txt'];

big_delta   = 25.5 * 1e-3;    %ms
small_delta = 10.8 * 1e-3;    %ms




%% Load data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Load subject information
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if verbose
%     disp('Loading config file...');
% end
% [ INFO ] = read_cfg_file( fconfig );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Load acquisition information
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if verbose
%     disp('Loading acq specs file...');
% end
% [ ACQ ] = read_acqspecs_file( fspecs );


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
DATAIN.lowb0 = DATAIN.dwi(:,:,:,1);

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

%%%% MEAN LOWB VOLUME
DATAIN.lowb = mean(DATAIN.dwi(:,:,:,index_b0),4);

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


if 1 == 1
%     qvol = shelldens(DATAIN.bvals);
    [qvol, res] = shelldens_rj(DATAIN.bvals);
    
    % display correction factor
    figure
    plot(res.bval_uniq, res.qvol_samp, '-o');
    grid on;
    ylim([0, 1]);
    title('sampling density nonuniformity correction factor');
    xlabel('b-value (s/mm^2)');
    ylabel('normalized correction factor');
    
else
    
    %% calculate q-space sampling density nonuniformity correction factor (Fig.4b in [1])

    bval = DATAIN.bvals(:); % column vector
    disp('double check unique b values are:');
    disp(unique(bval)); % note a +/- 50 variance at b=10k
    bval_rnd = round(bval / 200) * 200; % round up bval

    bval_uniq = unique(bval_rnd); % unique bval
    count = zeros(size(bval_uniq)); % number of each bval
    for ii = 1 : length(bval_uniq)
        count(ii) = sum(bval == bval_uniq(ii));
    end

    qval_uniq = sqrt(bval_uniq / max(bval_uniq)); % nomalize unique qval, q~sqrt(b)
    qval_contour = (qval_uniq(1 : end-1) + qval_uniq(2 : end)) / 2; % middle contours
    qval_contour = [qval_contour; 2 * qval_uniq(end) - qval_contour(end)]; % extrapolate one outer contour 

    qvol_shell = diff(qval_contour .^ 3); % qspace volume associated with each shell
    qvol_shell = [qval_contour(1) .^ 3; qvol_shell]; % add in central sphere volume

    qvol_shell = qvol_shell / qvol_shell(1); % normalize central sphere volume to 1
    qvol_samp = qvol_shell ./ count; % qspace volume associated with a sample on a shell

    qvol = zeros(size(bval)); % qspace volume for each sample (C in Eq.4 in [1])
    for ii = 1 : length(bval_uniq)
        b = bval_uniq(ii);
        qvol(bval == b) = qvol_samp(ii);
    end

    % or can be simply computed using shelldens.m
    % qvol = shelldens(bval, 200);

    % display correction factor
    figure
    plot(bval_uniq, qvol_samp, '-o');
    grid on;
    ylim([0, 1]);
    title('sampling density nonuniformity correction factor');
    xlabel('b-value (s/mm^2)');
    ylabel('normalized correction factor');

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Prepare for GDSI recon
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%
%%%%% NOTE: 
%%%%%   Orginally not sure what Delta/delta were, so used values from Fan
%%%%%   et al. (i think, or related bay 8 paper). 
%%%%%    - was using Delta/delta = 21.8/12.9 ms --> MDDwater = 16.20 um
%%%%%    - actually  Delta/delta = 25.5/10.8 ms --> MDDwater = 18.12 um    
%%%%%

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



% f=figure('color','w','position',[900 400 900 500],'InvertHardcopy','off'); 
% ax = axes('position',[0.05 0.05 0.4 0.9]);
% scatter3(nqvec(:,1),nqvec(:,2),nqvec(:,3), ...
%     30,vecnorm(nqvec,2,2),'filled'); 
% set(ax,'Colormap',jet);
% axis equal;
% ax = axes('position',[0.5 0.05 0.4 0.9]);
% scatter3(odf_verticesI(:,1),odf_verticesI(:,2),odf_verticesI(:,3), ...
%     30,vecnorm(odf_verticesI,2,2),'filled'); 
% set(ax,'Colormap',jet); caxis([0 1])
% axis equal;



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


%% For ODFs
if ReconODFs
    if verbose
        disp('Setting odf params...');
    end
    odfpar.idx_start = 1;
    odfpar.idx_end   = 60;
    odfpar.r_weight  = 2;
    odfpar.Wmtx      = repmat(dispvec.^odfpar.r_weight, ...
        nb_pdf_dirs, 1);
    odfpar.sep       = 25; %degrees, peak separation
    odfpar.clip      = true; %clip negative eap
    
    TR2 = triangulation((odf_faces+1)',odf_vertices');
end


%% For ZDPs
if ReconZDPs
    if verbose
        disp('Setting RTPP params...');
    end
    v = [1;0;0];
    a_circle = (0:pi/50:2*pi);
    a_circle = a_circle(2:end);
    a_plane = cat(1,zeros(size(a_circle)),cos(a_circle),sin(a_circle));
    % rtap_matrix = permute(repmat(rs(2:end),1,3,100),[3 2 1]);
    % [ rtap_Fmtx ] = precompute_Fmtx( qvec, a_plane, rs );
    zpds.v = v;
    zpds.a_plane = a_plane;
    zpds.qvec = qvec;
    zpds.rs = rs;
end
        

%% Precompute Fmtx
if verbose
    disp('Creating Fmtx...');
end
[ Fmtx ] = precompute_Fmtx( qvec, pdf_dirs, rs );          


%% Create output structure
if verbose
    disp('Creating output matrices...');
end
voldims = size(DATAIN.mask,1:2);
% [ DATAOUT, ~ ] = allocate_gdsi_v4( voldims, nr ); 
% if StoreRmatrix
%     DATAOUT.Rmatrix = zeros([voldims nr nb_pdf_dirs]);
% end

DATAOUT = [];

DATAOUT.Rmatrix = zeros([voldims nr nb_pdf_dirs]);

DATAOUT.Rvec = zeros([voldims nr]);
DATAOUT.Rstd = zeros([voldims nr]);
DATAOUT.Rrms = zeros([voldims nr]);
DATAOUT.MeanGA = zeros(voldims);
DATAOUT.MeanEAP = zeros(voldims);
DATAOUT.StdEAP = zeros(voldims);
DATAOUT.RmsEAP = zeros(voldims);
DATAOUT.EntropyEAP = zeros(voldims);

DATAOUT.MeanKurt = zeros(voldims);
DATAOUT.NG = zeros(voldims);
DATAOUT.KLDg = zeros(voldims);
DATAOUT.KLDi = zeros(voldims);
DATAOUT.NLrms = zeros(voldims);

DATAOUT.RTOP = zeros(voldims);
% DATAOUT.RTAP = zeros(voldims);
% DATAOUT.RTPP = zeros(voldims);
% DATAOUT.V1 = zeros([voldims 3]);

DATAOUT.HWHM = zeros(voldims);
DATAOUT.HWHMstd = zeros(voldims);
DATAOUT.HWHMmin = zeros(voldims);
DATAOUT.HWHMmax = zeros(voldims);
DATAOUT.HWHMmin_vec = zeros([voldims 3]);
DATAOUT.HWHMmax_vec = zeros([voldims 3]);

DATAOUT.qiv = zeros(voldims);
DATAOUT.qvar = zeros(voldims);
DATAOUT.qmsd = zeros(voldims);
DATAOUT.qrtop = zeros(voldims);
DATAOUT.qentropy = zeros(voldims);


%% BEGIN RECON!
if verbose
    disp(' --- Starting GDSI recon --- ');
end

%%%% % % % CHANGE MASK HERE!!!! % % % %%%%
DATAIN.mask(:,:,[1:34,36:end]) = 0;

totalvoxels = nnz(DATAIN.mask);

%%% Voxels of interest:
% 
% CS_DSI_10_01_2021: cc --- ix=65, iy=49, iz=35
% CS_DSI_10_01_2021: cross --- ix=56, iy=45, iz=35

cnt = 0;
ftic = tic;
rectimer = tic; timerstepsize = 100;
for ix=1:size(DATAIN.mask,1)
    for iy=1:size(DATAIN.mask,2)
        for iz=35 %1:size(DATAIN.mask,3)
            if DATAIN.mask(ix,iy,iz)>0
                
                if verbose && cnt>0 && mod(cnt,timerstepsize)==0
                    currtimer = toc(rectimer);
                    fprintf(' - %d/%d vox = %.02f pct\n',cnt,totalvoxels,(cnt/totalvoxels)*100);
                    fprintf('  - Elapsed time:             %.01f sec \n',currtimer);
                    secremaining = round((currtimer/cnt)*((totalvoxels-cnt)));
                    fprintf('  - Estimated time remaining: %d h %d m %d s \n',...
                        floor(secremaining/3600),...
                        floor(mod(secremaining,3600)/60),...
                        mod(secremaining,60));
                end

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % RECONSTRUCT EAP
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                % Get voxel qspace data as col vector
                data = reshape(DATAIN.dwi(ix,iy,iz,:),[],1);
                b0s = data(index_b0);       % b0 volumes
                dwis = data(~index_b0);     % dwi volumes
                data = [mean(b0s); dwis];   % average b0 volumes
                data = data / mean(b0s);    % normalize diff signal
                
                % metrics on q-space signal 
                DATAOUT.qiv(ix,iy) = inv(sum(data));
                DATAOUT.qmsd(ix,iy) = mean(data.*(mag_qvec.^2));
                DATAOUT.qvar(ix,iy) = var(data);
                DATAOUT.qrtop(ix,iy) = sum(data);
                DATAOUT.qentropy(ix,iy) = entropy(data);

                %%% Compute 1D EAP along each pdf_dir
                [ Rmatrix ] = compute_1d_eap_mtx( diag(qvol) * data, Fmtx );
                if StoreRmatrix
                    DATAOUT.Rmatrix(ix,iy,:,:) = Rmatrix;
                end
                
                Rvec        = mean(Rmatrix,2);
                Rstd        = std(Rmatrix,0,2);
                Rrms        = rms(Rmatrix,2);
                Rsigma      = std(Rmatrix(:));   %std of eap
                
                
                DATAOUT.RTOP(ix,iy)             = Rvec(1); %Rmatrix(1,1);      %rtop of eap
                
                DATAOUT.Rvec(ix,iy,:)           = Rvec;      
                DATAOUT.Rstd(ix,iy,:)           = Rstd;       
                DATAOUT.Rrms(ix,iy,:)           = Rrms;
                
                DATAOUT.MeanGA(ix,iy)           = mean(Rstd); 
                DATAOUT.MeanEAP(ix,iy)          = mean(Rmatrix(:));  %mean of eap
                DATAOUT.StdEAP(ix,iy)           = Rsigma;
                DATAOUT.RmsEAP(ix,iy)           = rms(Rmatrix(:));
                DATAOUT.EntropyEAP(ix,iy)       = entropy(Rmatrix(:));
                


                
                %%%%%%%% (NON)GAUSSIAN METRICS  %%%%%%%%
                % Divergence of Mean 1D EAP from Gaussian
                gaussmtx = exp(-dispvec2'*Rsigma);        %make gaussian
                gaussmtx = (gaussmtx/sum(gaussmtx))*sum(Rvec(:));
                                
                if min(Rvec)<0
                    rtmp = Rvec-min(Rvec);
                else
                    rtmp = Rvec;
                end
%                 if min(Rmatrix)<0
%                     mtmp = Rmatrix-min(Rmatrix(:));
%                 else
%                     mtmp = Rmatrix;
%                 end
                
                linearmtx = linspace(max(rtmp),0,nr)';
                
                isomtx = mean(rtmp)*ones(size(rtmp));
%                 isomtx = repmat(rtmp,1,size(Rmatrix,2));
                
%                 (KLDiv(mtmp(1:50,:),isomtx(1:50,:)));   %compute KL divergence
               
                DATAOUT.KLDg(ix,iy)   = (KLDiv(rtmp(1:50)',gaussmtx(1:50)'));   %compute KL divergence
                DATAOUT.KLDi(ix,iy)   = (KLDiv(rtmp(1:50)',isomtx(1:50)'));   %compute KL divergence
                DATAOUT.NLrms(ix,iy) = sqrt(mean((rtmp-linearmtx).^2));
                
                % Calculate non-Gaussianity (NG) [Ning 2015]
                costheta = dot(Rvec,gaussmtx)/sqrt(dot(Rvec,Rvec)*dot(gaussmtx,gaussmtx));
                DATAOUT.NG(ix,iy) = sqrt(1-costheta^2); %really sintheta
                
                %%%%%%%% KURTOSIS %%%%%%%%
                % Calculate kurtosis at each distance r
%                 DATAOUT.RotKurt(ix,iy,:) = kurtosis(Rmatrix');
                
                % Calculate average eap kurtosis
                DATAOUT.MeanKurt(ix,iy) = kurtosis(Rmatrix(:));
                
%                 % Mean squared displacement
%                 DATAOUT.MSD(ix,iy) = sum(dispvec2' .* Rvec);
%                 DATAOUT.MFD(ix,iy) = sum(dispvec4' .* Rvec);
                
                % HWHM radius (only half width, not full!)
                [ hwhm ] = compute_ralpha( Rmatrix, 0.5, dispvec );
                DATAOUT.HWHM(ix,iy) = hwhm.mean;
                DATAOUT.HWHMstd(ix,iy) = hwhm.std;
                DATAOUT.HWHMmin(ix,iy) = hwhm.min;
                DATAOUT.HWHMmax(ix,iy) = hwhm.max;
                
                tmp_mvec = pdf_dirs((hwhm.vec==hwhm.min),:);
                if size(tmp_mvec,1)>1
                    tmp_mvec = tmp_mvec(1,:);
                end
                DATAOUT.HWHMmin_vec(ix,iy,:) = tmp_mvec;
                
                tmp_mvec = pdf_dirs((hwhm.vec==hwhm.max),:);
                if size(tmp_mvec,1)>1
                    tmp_mvec = tmp_mvec(1,:);
                end
                DATAOUT.HWHMmax_vec(ix,iy,:) = tmp_mvec;
                
                if ClipRmatrix
                    % Clip after first zero crossing
                    [ nRmatrix ] = clip_negative_Rmatrix( Rmatrix );
%                     [ rRmatrix ] = clip_ringing_Rmatrix( Rmatrix );
%                     [ nrRmatrix ] = clip_negative_Rmatrix( rRmatrix );
                    nRvec = mean(nRmatrix, 2);
                    nRstd = std( nRmatrix, 0, 2);

                    if MakePlots.rgbprofiles
                        [ f1 ] = plot_all_profiles_voxel( dispvec, Rmatrix, ...
                            abs(diff_disp_sphere.vertices)' );
                         hold on; plot(dispvec,zeros(size(dispvec)),'-k');
                         
                        pause(1);
                        
                        [ f1 ] = plot_all_profiles_voxel( dispvec, nRmatrix, ...
                            abs(diff_disp_sphere.vertices)' ); 
                         hold on; plot(dispvec,zeros(size(dispvec)),'-k');
                    end
                end
                

                
                if ReconODFs
                    %%%%%%%% ODF %%%%%%%%
                    % Recon dODF
                    Omatrix = Rmatrix.*odfpar.Wmtx';
                    if Use_DSI_Studio_dirs
                        Omatrix = repmat(Omatrix,1,2);
                    end                
                    ODF = sum( Omatrix( odfpar.idx_start:odfpar.idx_end, :) );
                    ODF_points = diff_disp_sphere.vertices .* repmat(ODF', [1, 3]);
                    
                    [ f1 ] = plot_all_profiles_voxel( dispvec, Omatrix, ...
                            abs(diff_disp_sphere.vertices)' ); 
                         hold on; plot(dispvec,zeros(size(dispvec)),'-k');
                         
                         
                         
                    nOmatrix = nRmatrix.*odfpar.Wmtx';
                    if Use_DSI_Studio_dirs
                        nOmatrix = repmat(nOmatrix,1,2);
                    end                
                    nODF = sum( nOmatrix ); %( odfpar.idx_start:odfpar.idx_end, :) );
                    nODF_points = diff_disp_sphere.vertices .* repmat(nODF', [1, 3]);
                    
                    [ f1 ] = plot_all_profiles_voxel( dispvec, nOmatrix, ...
                            abs(diff_disp_sphere.vertices)' ); 
                         hold on; plot(dispvec,zeros(size(dispvec)),'-k');
                         
                    [pks2, allpks2] = gdsi_find_peak(nODF', diff_disp_sphere, odfpar.sep);     
                         

                    % find ODF peaks
                    [pks, allpks] = gdsi_find_peak(ODF', diff_disp_sphere, odfpar.sep);
                    npeaks.D19 = length(pks.D19.pkvals);
                    totalpeaks.D19 = npeaks.D19;
                    if npeaks.D19>3                    
                        npeaks.D19=3;
                    end
                    
                    %%% Compute RTAP, RTPP on major fiber axis
                    pdd = pks.pkvecs(1,:);
                    [ rtap ] = compute_RTAP( data, v, pdd, qvec, rs, dispvec2 );
                    [ rtpp ] = compute_RTPP( data, v, a_plane, pdd, qvec, rs, dispvec2 );

                    %%% Compute fixel metrics
                    [fixels] = compute_fixel_metrics_v2( pks, diff_disp_sphere.vertices, ...
                        Rmatrix, data, dispvec, zpds); 


                    if MakePlots.odf
                        %%%%% plot ODF
    %                     plotODF = ODF;
                        plotODF = ODF/max(ODF);
        % %                 plotODF = (ODF-min(ODF))/(max(ODF)-min(ODF));
                        plotPeaks = pks.pkvecs;

                        [ f ] = gdsi_plot_odf( diff_disp_sphere, plotODF, plotPeaks );
                        
                        
                        
                        
                        plotODF = nODF/max(nODF);
        % %                 plotODF = (ODF-min(ODF))/(max(ODF)-min(ODF));
                        plotPeaks = pks2.pkvecs;

                        [ f ] = gdsi_plot_odf( diff_disp_sphere, plotODF, plotPeaks );
                        
                        
                    end
                
                    %%% STORE RESULTS
                    % ODF GFA, RMS, STD, npeaks, iso, mass, max
                    DATAODF.GFA(ix,iy) = std(nODFhalf)/rms(nODFhalf);
                    DATAODF.RMS(ix,iy) = rms(nODFhalf);
                    DATAODF.STD(ix,iy) = std(nODFhalf);
                    DATAODF.npeaks(ix,iy) = n_totalpeaks;
                    % ODF iso (min value)
                    DATAODF.iso(ix,iy) = min(nODF);
                    % ODF mass (sum of all odf values)
                    DATAODF.mass(ix,iy) = sum(nODF);

                    %%% Fixels
                    for nn=1:n_npeaks
                        eval(['curr = n_fixels.v' num2str(nn) ';']);
                        currind = curr.diridx;
                        currO = nOmatrix(odfpar.idx_start:odfpar.idx_end, currind);
                        currfa = sum(currO);
    %                     currnfa = currfa/sum(currO);

                        eval(['DATAODF.ind' num2str(nn) '(ix,iy) = curr.diridx;']);
                        eval(['DATAODF.vec' num2str(nn) '(ix,iy,:) = curr.dir;']);

                        eval(['DATAODF.f' num2str(nn) '(ix,iy) = curr.f;']);
                        eval(['DATAODF.fa' num2str(nn) '(ix,iy) = curr.fa;']);

                        eval(['DATAODF.Rvec' num2str(nn) '(ix,iy,:) = curr.Rvec;']);
                        eval(['DATAODF.rtap' num2str(nn) '(ix,iy) = curr.Rmean;']);
                        eval(['DATAODF.Rstd' num2str(nn) '(ix,iy) = curr.Rstd;']);
                        eval(['DATAODF.Rmsd' num2str(nn) '(ix,iy) = curr.Rmsd;']);            
                        eval(['DATAODF.kurt' num2str(nn) '(ix,iy) = curr.Kurtosis;']);                    
                        eval(['DATAODF.ip' num2str(nn) '(ix,iy) = curr.inflpnt;']);
                        eval(['DATAODF.r_50_' num2str(nn) '(ix,iy) = curr.r_50;']);
                        eval(['DATAODF.r_25_' num2str(nn) '(ix,iy) = curr.r_25;']);

                        eval(['DATAODF.rtpp.Rvec' num2str(nn) '(ix,iy,:) = curr.rtpp.Rvec;']);                    
                        eval(['DATAODF.rtpp.Rstd' num2str(nn) '(ix,iy,:) = curr.rtpp.Rstd;']);                    
                        eval(['DATAODF.rtpp.rtpp' num2str(nn) '(ix,iy,:) = curr.rtpp.RTPP;']);                    
                        eval(['DATAODF.rtpp.rtpv' num2str(nn) '(ix,iy,:) = curr.rtpp.RTPV;']);                    
                        eval(['DATAODF.rtpp.msd' num2str(nn) '(ix,iy,:) = curr.rtpp.msd;']);                    
                        eval(['DATAODF.rtpp.kurt' num2str(nn) '(ix,iy,:) = curr.rtpp.kurtosis;']);                    

                    end
                end
                

                
                cnt=cnt+1;
            end
        end
    end
end

% foutdir = '/autofs/space/hemera_002/users/rjjones/gdap/subjects/qt';
% foutmat = [foutdir filesep 'z0037_Rmatrix.mat'];
% save(foutmat,'DATAOUT');


disp(' --- Finished reconstructing voxels --- ');
ftoc=toc(ftic);
fprintf('%.1f seconds | %d voxels\n',ftoc,totalvoxels);


disp(' - Saving output .mat file...');
% outdir = [procdir filesep 'gdsi_recon_vol'];

if ~exist(outdir,'dir')
    fprintf('Making outdir:\n %s\n',outdir);
    mkdir(outdir);
end

niioutdir = [outdir filesep 'gdsi_recon_csdsi-nufft-hcp3sh5k_with-shelldens-corr_021022_slice_z35'];
if ~exist(niioutdir,'dir')
    mkdir(niioutdir);
end



bvals               = DATAIN.bvals;
bvecs               = DATAIN.bvecs;
nodif_brain_mask    = DATAIN.mask;
lowb                = DATAIN.lowb;

% foutmat = [outdir filesep 'gdsi_recon_csdsi_subcortical-only.mat'];
foutmat = [outdir filesep 'gdsi_recon_csdsi-nufft-hcp3sh5k_with-shelldens-corr_021022_slice_z35.mat'];

fprintf('Saving recon to file:\n %s\n',foutmat);
save(foutmat, 'DATAOUT', ... %'INFO', 'ACQ',...
   'pdf_dirs', 'dispvec', ...
    'diff_disp_sphere', 'qvec', 'D_water', 'mdd_water', ...
    'bvals', 'bvecs', 'totalvoxels', 'ftoc', ...
    'nodif_brain_mask', '-v7.3');


fprintf('\n-------------\n - Finished GDSI recon - \n-------------\n');



dec_hwhm_min = abs(DATAOUT.HWHMmin_vec).*repmat((DATAOUT.HWHMmin-4)/1.5,1,1,3);
figure; imshow(dec_hwhm_min); title('DEC HWHM_{min}','Interpreter','tex');

dec_hwhm_max = abs(DATAOUT.HWHMmax_vec).*repmat((DATAOUT.HWHMmax-5.5)/2,1,1,3);
figure; imshow(dec_hwhm_max); title('DEC HWHM_{max}','Interpreter','tex');


f=figure('color','w','InvertHardcopy','off','position',[500 800 400 400]);
imshow(dec_hwhm_min,'border','tight');
fout = [niioutdir filesep 'DEC_HWHM_min_z35.png'];
print(f,fout,'-dpng');
clf;
imshow(dec_hwhm_max,'border','tight');
fout = [niioutdir filesep 'DEC_HWHM_max_z35.png'];
print(f,fout,'-dpng');
clf;



%%%%%%%%%%%%

vol=load_untouch_nii('/autofs/space/nyx_002/users/rjones/Bay8/CS_DSI_10_01_2021/preproc/anat2diff/fs-tissue-masks/all-wmmask.nii.gz');
wmmask=double(vol.img);
wmslice = (wmmask(:,:,35));

vol=load_untouch_nii('/autofs/space/nyx_002/users/rjones/Bay8/CS_DSI_10_01_2021/preproc/anat2diff/fs-tissue-masks/ctx-gmmask.nii.gz');
gmmask=double(vol.img);
gmslice = (gmmask(:,:,35));

curr = DATAOUT.Rvec;
for rr=1:nr
    tmp=curr(:,:,rr);
    tmp(wmslice==0)=0;
    a=mean(tmp(tmp>0));
    b=std(tmp(tmp>0));
    dat.rvec.wm.mean(rr)=a;
    dat.rvec.wm.std(rr)=b;
end

curr = DATAOUT.Rvec;
for rr=1:nr
    tmp=curr(:,:,rr);
    tmp(gmslice==0)=0;
    a=mean(tmp(tmp>0));
    b=std(tmp(tmp>0));
    dat.rvec.gm.mean(rr)=a;
    dat.rvec.gm.std(rr)=b;
end

curr = DATAOUT.Rstd;
for rr=1:nr
    tmp=curr(:,:,rr);
    tmp(wmslice==0)=0;
    a=mean(tmp(tmp>0));
    b=std(tmp(tmp>0));
    dat.rstd.wm.mean(rr)=a;
    dat.rstd.wm.std(rr)=b;
end

curr = DATAOUT.Rstd;
for rr=1:nr
    tmp=curr(:,:,rr);
    tmp(gmslice==0)=0;
    a=mean(tmp(tmp>0));
    b=std(tmp(tmp>0));
    dat.rstd.gm.mean(rr)=a;
    dat.rstd.gm.std(rr)=b;
end

figure; hold on; plot(dat.rvec.wm.mean); plot(dat.rvec.gm.mean);title('Mean EAP - hcp 3sh 5k');
figure; hold on; plot(dat.rstd.wm.mean); plot(dat.rstd.gm.mean);title('GA profile - hcp 3sh 5k');






% % Save volumes in .nii format
fafile='/autofs/space/nyx_002/users/rjones/Bay8/CS_DSI_10_01_2021/preproc/dtifit/dtifit_FA.nii.gz';
vol = load_untouch_nii(fafile);
% vol.hdr.dime.datatype = 16;
% vol.hdr.dime.bitpix = 32;





filelist = {'RTOP','StdEAP','MeanEAP','MeanGA','MeanKurt','KLDg',...
            'KLDi','NLrms','EntropyEAP','RmsEAP','NG',...
            'HWHM','HWHMstd','HWHMmin','HWHMmax',...
            'qiv','qvar','qmsd','qentropy'};
for f=1:length(filelist)
    F = filelist{f};
    eval(['curr=DATAOUT.' F ';']);
    tmp = zeros(size(vol.img));
    tmp(:,:,35)=curr;
    outfile = [niioutdir filesep  F '.nii.gz'];
    write2nifti( vol, tmp, outfile );
end

filelist = {'HWHMmin_vec','HWHMmax_vec'  }; 
for f=1:length(filelist)
    F = filelist{f};
    eval(['curr=DATAOUT.' F ';']);
    tmp = zeros([size(vol.img) 3]);
    tmp(:,:,35,:)=curr;
    outfile = [niioutdir filesep  F '.nii.gz'];
    write2nifti( vol, tmp, outfile );
end


filelist = {'Rvec','Rstd'};
for f=1:length(filelist)
    F = filelist{f};
    eval(['curr=DATAOUT.' F ';']);
    tmp = zeros([size(vol.img) 100]);
    tmp(:,:,35,:) = curr;
    outfile = [niioutdir filesep  F '.nii.gz'];
    write2nifti( vol, tmp, outfile );
end

