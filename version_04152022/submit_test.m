% [ DATAOUT ] = script_GDSI_recon_test( fdir, odir, isshelled, verbose ) 

dtroot = '/autofs/space/nyx_002/users/rjones/Bay8/CS_DSI_10_01_2021/preproc/nufft';
verbose = 1;

fdir = [dtroot filesep '3sh'];
odir = [fdir filesep 'gdsi'];
isshelled = true;
script_GDSI_recon_test( fdir, odir, isshelled, verbose );

% fdir = [dtroot filesep '3sh'];
% odir = [fdir filesep 'gdsi'];
% isshelled = true;
% script_GDSI_recon_test( fdir, odir, isshelled, verbose );
% 
% fdir = [dtroot filesep '2sh'];
% odir = [fdir filesep 'gdsi'];
% isshelled = true;
% script_GDSI_recon_test( fdir, odir, isshelled, verbose );
% 
% fdir = [dtroot filesep 'FS-DSI'];
% odir = [fdir filesep 'gdsi'];
% isshelled = false;
% script_GDSI_recon_test( fdir, odir, isshelled, verbose );



