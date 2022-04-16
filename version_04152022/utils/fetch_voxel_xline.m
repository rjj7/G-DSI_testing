function [ xline ] = fetch_voxel_xline( invol, y, z )
% [ xline ] = fetch_voxel_xline( invol, y, z )

xline = squeeze(invol(:,y,z));
