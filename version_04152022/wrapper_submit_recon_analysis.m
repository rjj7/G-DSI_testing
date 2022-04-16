
dtroot      = '/autofs/space/nyx_002/users/rjones/Bay8/CS_DSI_10_01_2021';

filedir     = [dtroot filesep 'diff/preproc/mri/eddy/results_def/nufft/hcp_2sh_1k,3k'];
GDSI_recon_eddy_def_nufft( dtroot, filedir );
script_GDSI_analysis_nufft( dtroot, filedir );

filedir     = [dtroot filesep 'diff/preproc/mri/eddy/results_def/nufft/hcp_2sh_1k,5k'];
GDSI_recon_eddy_def_nufft( dtroot, filedir );
script_GDSI_analysis_nufft( dtroot, filedir );

filedir     = [dtroot filesep 'diff/preproc/mri/eddy/results_def/nufft/hcp_2sh_1k,7k'];
GDSI_recon_eddy_def_nufft( dtroot, filedir );
script_GDSI_analysis_nufft( dtroot, filedir );

