function [ cRmatrix ] = clip_ringing_Rmatrix( Rmatrix )
    
% clip ringing (after any increase in 1d eap profile)
cRmatrix = Rmatrix;
a = diff(cRmatrix);
for dd=1:size(a,2)
    cur = cRmatrix(:,dd);
    curr = a(:,dd);
    idx=find(curr>0,1,'first');
    if ~isempty(idx) && idx<length(curr)
        cur(idx+1:end) = 0;
    end
    cRmatrix(:,dd) = cur;
end