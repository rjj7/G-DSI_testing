function [ DATAOUT ] = allocate_gdsi( voldims, nr, ndirs, nr_odf )
%
% [ DATAOUT ] = allocate_gdsi( voldims, nr );
% 
%  DATAOUT.sph_avg_mean
%  DATAOUT.sph_avg_mean_nneg
%  DATAOUT.sph_avg_stdev
%  
% 

if nargin<2
    disp('USAGE: [ DATAOUT ] = allocate_gdsi( voldims, nr )');
    return
end
if nargin<3
    ndirs=[];
end
if nargin<4
    nr_odf=[];
end

DATAOUT = [];

%Spherically averaged prob. dens. at each displ. radius
DATAOUT.avg_disp_sph        = zeros([voldims nr]);
% DATAOUT.avg_disp_sph_nneg   = zeros([voldims nr]);
DATAOUT.std_disp_sph        = zeros([voldims nr]);

%For KL divergence vs. isotropic
DATAOUT.isokldiv            = zeros([voldims nr]);
DATAOUT.eap_kurtosis        = zeros([voldims nr]);
if ~isempty(ndirs)
    DATAOUT.rad_kurtosis        = zeros([voldims ndirs]);
else
    DATAOUT.rad_kurtosis = [];
end
if ~isempty(nr_odf)
    DATAOUT.odfs = zeros([voldims nr_odf]);
    DATAOUT.v1 = zeros([voldims 3]);
    DATAOUT.v2 = zeros([voldims 3]);
    DATAOUT.v3 = zeros([voldims 3]);
    DATAOUT.odf_iso = zeros(voldims);
    DATAOUT.f1 = zeros(voldims);
    DATAOUT.f2 = zeros(voldims);
    DATAOUT.f3 = zeros(voldims);
    DATAOUT.msd1 = zeros(voldims);
    DATAOUT.msd2 = zeros(voldims);
    DATAOUT.msd3 = zeros(voldims);
    
end


% %Std. Dev. across sphere/dirs of prob. dens. at each displ. radius
% DATAOUT.var_disp_sph        = zeros([voldims nr]);
% DATAOUT.var_disp_sph_nneg   = zeros([voldims nr]);

%Spherically averaged R_0 at each voxel
% DATAOUT.r0_indices_sph      = zeros(voldims);

% %R_0 at each dir in each voxel
% DATAOUT.r0_indices_ALL      = zeros([voldims ndir]);

end