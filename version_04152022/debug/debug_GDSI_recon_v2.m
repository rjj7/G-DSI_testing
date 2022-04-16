function [ DATAOUT, res ] = debug_GDSI_recon_v2( DATAIN, xrange, yrange, zrange, Fmtx ) 
% 
% Time GDSI recon, using concatenated Fmtx (DR x N)
%
% original script by (c) Qiyuan Tian, Stanford RSL
% Modified by RJ, Feb 2021 - Feb 2022
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isdeployed
    disp('--RUNNING DEPLOYED APP--');
else
    addpath('/autofs/space/hemera_001/users/rjjones/code/toolboxes/NIfTI_20140122');
    addpath(genpath('/autofs/space/nyx_002/users/rjones/Bay8/code/gdsi/version_020122'));
end

% clear; close all; clc;

DATAOUT = [];
% DATAOUT.Rmatrix = zeros([voldims nr nb_pdf_dirs]);

res = [];

disp(' --- Starting GDSI recon --- ');

cnt = 0;
ftic = tic;
for inx=1:length(xrange)
    ix=xrange(inx);
    for iny=1:length(yrange )
        iy=yrange(iny);
        for inz=1:length(zrange)
            iz=zrange(inz);
            if DATAIN.mask(ix,iy,iz)>0
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % RECONSTRUCT EAP
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                % Get voxel qspace data as col vector
                data = reshape(DATAIN.dwi(ix,iy,iz,:),[],1);
                b0s = data(DATAIN.index_b0);       % b0 volumes
                dwis = data(~DATAIN.index_b0);     % dwi volumes
                data = [mean(b0s); dwis];   % average b0 volumes
                data = data / mean(b0s);    % normalize diff signal
                
%                 %%% Compute 1D EAP along each pdf_dir
%                 [ Rmatrix ] = compute_1d_eap_mtx( data, Fmtx );
                
                Rmatrix = Fmtx * data;

%                 if StoreRmatrix
%                     DATAOUT.Rmatrix(ix,iy,:,:) = Rmatrix;
%                 end

                cnt=cnt+1;
            end
        end
    end
end
ftoc=toc(ftic);

fprintf('%.1f seconds | %d voxels\n',ftoc,cnt);

res.ftoc = ftoc;
res.cnt = cnt;

% disp(' - Saving output .mat file...');
% % outdir = [procdir filesep 'gdsi_recon_vol'];
% 
% % foutmat = [outdir filesep 'gdsi_recon_csdsi_subcortical-only.mat'];
% foutmat = [outdir filesep 'gdsi_recon_csdsi_021022_slice_z35.mat'];
% 
% fprintf('Saving recon to file:\n %s\n',foutmat);
% save(foutmat, 'res', '-v7.3');


% fprintf('\n-------------\n - Finished GDSI recon - \n-------------\n');
