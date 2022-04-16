function [ nRmatrix ] = clip_negative_Rmatrix( Rmatrix )

nRmatrix = Rmatrix;
nRmatrix(Rmatrix<=0) = 0;
[~,b]=find(nRmatrix==0);
d=unique(b);
for dd=1:length(d)
    d_=d(dd);
    cur=nRmatrix(:,d_);
    idx=find(cur==0,1,'first');
    cur(idx+1:end)=0;
    nRmatrix(:,d_)=cur;
end