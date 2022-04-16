function [ ralpha ] = compute_ralpha( Rmatrix, p, dispvec )
% [ ralpha ] = compute_ralpha( Rmatrix, p, dispvec )
% 
%  Rmatrix is stack of 1D EAP profiles
%  p is target prob. dens. to use
%  dispvec is displacement distance vector
% 
% ralpha is struct with info on radius where EAP(ralpha)=p

nd=size(Rmatrix,2);
ralpha_vec=zeros(nd,1);
% r_step_size=diff(dispvec);
r_step_size=dispvec(2)-dispvec(1);
for dd=1:nd
    cur=Rmatrix(:,dd);
    cur=cur/cur(1);
    idx=find(cur<=p,1,'first');
    if ~isempty(idx)
        x1=dispvec(idx-1);
        y1=cur(idx-1);
        y2=cur(idx);
        deltax=r_step_size*((y1-p)/(y1-y2));
        xhat=x1+deltax;
        ralpha_vec(dd)=xhat;
    else
        ralpha_vec(dd)=NaN;
    end
end
% figure; plot(ralpha_vec); title(['r_{' num2str(ralpha) '}']);
% ralpha_vec(isnan(ralpha_vec))=[];
ralpha.vec  = ralpha_vec;
ralpha.vec(isnan(ralpha.vec))=0;
ralpha.mean = mean(ralpha_vec(~isnan(ralpha_vec)));
ralpha.std  = std(ralpha_vec(~isnan(ralpha_vec)));
ralpha.min = min(ralpha_vec(~isnan(ralpha_vec)));
ralpha.max = max(ralpha_vec(~isnan(ralpha_vec)));
ralpha.rng = ralpha.max-ralpha.min;