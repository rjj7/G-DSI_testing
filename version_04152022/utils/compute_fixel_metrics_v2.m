function [fixels] = compute_fixel_metrics_v2( pks, odf_vertices, Rmatrix, data, dispvec, zpds )

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
    
    currRvec   = Rmatrix(:,idx);
    currRvecSq = currRvec.*sq_dispvec;
%     currOvec   = currRvecSq;
%     currOvec   = Omatrix(:,idx);
%     currRvecFourth = currRvec.*fourth_dispvec;
    
    %%% Compute RTAP, RTPP
    [ rtpp ] = compute_RTPP( data, zpds.v, zpds.a_plane, curr, zpds.qvec, zpds.rs, sq_dispvec );
%     [ rtpp ] = compute_RTPP( data, zpds.v, curr, zpds.qvec, zpds.rs, sq_dispvec );
    
    
    %%%%%%%%% STORE RESULTS %%%%%%%%%
    %%% Direction/vector
    eval(['fixels.v' num2str(n) '.diridx = find(idx);']);
    eval(['fixels.v' num2str(n) '.dir = curr;']);
    
    %%% Vol. fraction [f = norm so f1=1; fa = raw odf values]
    eval(['fixels.v' num2str(n) '.f = pks.pkfracs(n);']);
    eval(['fixels.v' num2str(n) '.fa = pks.pkvals(n);']);
    
    %%% Fixel 1D EAP profile
    eval(['fixels.v' num2str(n) '.Rvec  = currRvec;']);
    
    %%% Fixel mean of 1D EAP
    eval(['fixels.v' num2str(n) '.Rmean = mean(currRvec);']);
    
    %%% Fixel std of 1D EAP 
    eval(['fixels.v' num2str(n) '.Rstd   = std(currRvec);']);
    
    %%% Fixel MSD
    eval(['fixels.v' num2str(n) '.Rmsd   = mean(currRvecSq);']);
    
    %%% Fixel std of squared 1D EAP 
    eval(['fixels.v' num2str(n) '.sqRstd = std(currRvecSq);']);
    
    %%% Fixel kurtosis of 1D EAP
    eval(['fixels.v' num2str(n) '.Kurtosis = kurtosis(currRvec);']);    
    
    %%% Fixel RTAP
    eval(['fixels.v' num2str(n) '.rtpp.Rvec = rtpp.Rvec;']);    
    eval(['fixels.v' num2str(n) '.rtpp.Rstd = rtpp.Rstd;']);
    eval(['fixels.v' num2str(n) '.rtpp.RTPP = rtpp.RTPP;']);
    eval(['fixels.v' num2str(n) '.rtpp.RTPV = rtpp.RTPV;']);
    eval(['fixels.v' num2str(n) '.rtpp.msd = rtpp.msd;']);
    eval(['fixels.v' num2str(n) '.rtpp.kurtosis = rtpp.kurtosis;']);
    
    
    %%% Fixel FWHM prob. dens.
    palphas = [0.5 0.25];
    for ra=1:length(palphas)
        p = palphas(ra);
        pstring = num2str(p*100);
        [ ralpha ] = compute_ralpha( currRvec, p, dispvec );
        eval(['fixels.v' num2str(n) '.r_' pstring ' = ralpha.mean;']);
        clear ralpha
    end

    
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
    eval(['fixels.v' num2str(n) '.inflpnt = xhat;']);
    clear x1 y1 y2 deltax
    
end

