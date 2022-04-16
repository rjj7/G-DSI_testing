function [ data_masked ] = mask_with_FS_label( data, labels, labelnb )
%
% [ data_masked ] = mask_with_FS_label( data, labels, labelnb )
%
% Mask input data with labels==labelnb
%  - data is input to mask
%  - labels is 3d vol of labels
%  - labelnb is vector of label nbs to use
% 

% find nb of vols of input data
nt = size(data,4);

% find nb of labels
nlabels = length(labelnb);

% get inds for label #1
islabeled = labels==labelnb(1);

% if more than 1 label, ...
if nlabels>1
    for l=2:nlabels
        tmp = labels==labelnb(l);
        islabeled = islabeled | tmp;
    end
end

% find total voxels
nvoxels = nnz(islabeled);

% get masked data
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
