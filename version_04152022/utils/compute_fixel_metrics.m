function [fixels] = compute_fixel_metrics( pks, odf_vertices, Rmatrix, Omatrix, dispvec )

fixels=[];

npks = length(pks.pkvals);
if npks<1
    disp(' NO PEAKS ');
    return
end

if size(odf_vertices,2)>size(odf_vertices,1)
    odf_vertices = odf_vertices';
end

if size(dispvec,1)<size(dispvec,2)
    dispvec=dispvec';
end

sq_dispvec = dispvec.^2;
% fourth_dispvec = dispvec.^4;

r_step_size = dispvec(2)-dispvec(1);

for n=1:npks
    curr = pks.pkvecs(n,:);
    idx = (odf_vertices(:,1)==curr(1) & ...
        odf_vertices(:,2)==curr(2) & ...
        odf_vertices(:,3)==curr(3));
    
    currRvec = Rmatrix(:,idx);
    currRvecSq = currRvec.*sq_dispvec;
    currOvec = Omatrix(:,idx);
%     currRvecFourth = currRvec.*fourth_dispvec;
    
    
%     [ rtap ] = compute_RTAP( data, v, a_plane, curr, qvec, rs, sq_dispvec );
%     [ rtpp ] = compute_RTPP( data, v, curr, qvec, rs, sq_dispvec );
                
                
    eval(['fixels.v' num2str(n) '.diridx = find(idx);']);
    eval(['fixels.v' num2str(n) '.dir = curr;']);
    eval(['fixels.v' num2str(n) '.f = pks.pkfracs(n);']);
    eval(['fixels.v' num2str(n) '.fa = pks.pkvals(n);']);
    
    eval(['fixels.v' num2str(n) '.msd = mean(currRvecSq);']);
    eval(['fixels.v' num2str(n) '.odf_msd = mean(currOvec);']);
    eval(['fixels.v' num2str(n) '.Rvec = currRvec;']);
%     eval(['fixels.v' num2str(n) '.sqRvec = currRvecSq;']);
%     eval(['fixels.v' num2str(n) '.Ovec = currOvec;']);
    
    eval(['fixels.v' num2str(n) '.Rstd = std(currRvec);']);
    eval(['fixels.v' num2str(n) '.Rsum = sum(currRvec);']);
    eval(['fixels.v' num2str(n) '.sqRstd = std(currRvecSq);']);
    eval(['fixels.v' num2str(n) '.sqRsum = sum(currRvecSq);']);
%     eval(['fixels.v' num2str(n) '.Ostd = std(currOvec(currOvec>0));']);
%     eval(['fixels.v' num2str(n) '.Osum = sum(currOvec(currOvec>0));']);
    
    % Fixel FWHM prob. dens.
    ind_fwhm = find(currRvec<(max(currRvec)/2),1,'first');
    %%% old way - not precise
%     rad_fwhm = dispvec(ind_fwhm);
    %%% new way - linear approx.
    x1=dispvec(ind_fwhm-1);
    y1=currRvec(ind_fwhm-1);
    y2=currRvec(ind_fwhm);
    deltax=r_step_size*((y1)/(y1-y2));
    rad_fwhm=x1+deltax;
    eval(['fixels.v' num2str(n) '.r_fwhm = rad_fwhm;']);
    clear x1 y1 y2 deltax
    
    % Fixel kurtosis
    eval(['fixels.v' num2str(n) '.MKPr = kurtosis(currRvec);']);
    
    
    %%%% Calculate inflection point
%     dr1=diff(fixels.v1.Rvec,1);
    dr2=diff(currRvec,2);

    if dr2(1)<0
        tmp=dr2>0;
        ind=find(tmp==1,1,'first');
    else
        tmp=dr2<0;
        ind=find(tmp==1,1,'first');
    end
    
    x1=dispvec(ind-1);
    y1=dr2(ind-1);
    y2=dr2(ind);
    deltax=r_step_size*((y1)/(y1-y2));
    xhat=x1+deltax;
%     xhat = dispvec(ind);
    
    eval(['fixels.v' num2str(n) '.inflpnt = xhat;']);
    clear x1 y1 y2 deltax
    
end

