function [ data_masked ] = mask_with_FS_trc( data, labels )
%
% [ data_masked ] = mask_with_FS_trc( data, labels )
%
% Mask input data with trc tractography logical mask
%  - labels is logical array for tract
% 

nt = size(data,4);

islabeled = labels==1;
nvoxels = nnz(islabeled);
data_masked = [];

if nvoxels>0
    data_masked = zeros(nvoxels,nt);
    for t=1:nt
        tmp = data(:,:,:,t);
        tmp = tmp(islabeled==1);
        data_masked(:,t) = tmp;
    end
else
    fprintf(' --- Could not find any voxels with label nb %d!!\n',labelnb);
end

end