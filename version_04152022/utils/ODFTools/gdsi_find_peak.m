function [pks, allpks] = gdsi_find_peak(odf, odf_patch, peaksep)
%     [pks, allpks] = gdsi_find_peak(odf, odf_patch[, peaksep])
%
%     INPUTS:  
%       odf = 1D ODF values
%       odf_patch = has fields [odf_patch.faces, odf_patch.vertices]
%       (opt.) peaksep = in degrees; default=25

% ~~~ From DSI Studio
% RJ , GDSI , 6/30/2021

if nargin<3
    peaksep = 25;
end

is_peak = odf;
maxodf = max(odf);
promthr = maxodf*1.05;

odf_faces = odf_patch.faces;
odf_faces = odf_faces + 1;

odf_vertices = odf_patch.vertices;

if size(odf_faces,1)>size(odf_faces,2)
    odf_faces = odf_faces';
end
if size(odf_vertices,1)>size(odf_vertices,2)
    odf_vertices = odf_vertices';
end

odf_faces = odf_faces - (odf_faces > length(odf))*length(odf);

is_peak(odf_faces(1,odf(odf_faces(2,:)) >= odf(odf_faces(1,:)) | ...
    odf(odf_faces(3,:)) >= odf(odf_faces(1,:)))) = 0;

is_peak(odf_faces(2,odf(odf_faces(1,:)) >= odf(odf_faces(2,:)) | ...
    odf(odf_faces(3,:)) >= odf(odf_faces(2,:)))) = 0;

is_peak(odf_faces(3,odf(odf_faces(2,:)) >= odf(odf_faces(3,:)) | ...
    odf(odf_faces(1,:)) >= odf(odf_faces(3,:)))) = 0;

[values,ordering] = sort(-is_peak);


% %% %  Old dsi studio:  % %% %
% p.inds = ordering(values < 0);
% p.vals = values(values<0);


% %% %  W/ peak sep thr:  % %% %
% peaks have negative values!
nvalues = find(values<0);
pks.sep = peaksep;
allpks = [];

% clear pks allpks
% If a peak was found, then:
if ~isempty(nvalues) 
%     n               = length(nvalues);          % # of peaks
    allpks.inds     = ordering(nvalues);        % inds of peaks
    allpks.inds = mod(allpks.inds,length(odf_vertices)/2);
    [uv,ui] = unique(allpks.inds,'stable');
    uv(uv==0)=length(odf_vertices)/2;
    allpks.inds = uv;
    n = length(allpks.inds);
    allpks.pks      = odf_vertices(:, ...
                         allpks.inds);          % all peaks detected
    allpks.vals     = -values(ui);         % values of peaks
    allpks.nvals    = allpks.vals/max(odf);            % norm. values of peaks
    allpks.seps     = zeros(size(allpks.pks));  % to store peak seps
    allpks.used     = zeros(n,1);               % store whether pk was used
    allpks.used(1)  = 1;
    
    pks.pkvecs      = allpks.pks(:,1);          % save first peak
    pks.pkfracs     = allpks.nvals(1);          % save first peak norm. val
    pks.pkvals      = allpks.vals(1);           % save first peak val
    % if g.t. 1 peak found
    if n>1 && size(allpks.pks,2)>1        
        % iterate over peaks 2:total
        for m=2:n
            % get current peak vector
            curr = allpks.pks(:,m);            
            % inds to all previous peaks
            h = 1:m-1;             
            % check sep vs all previous peaks
            sepflag = 1;
            for g=1:length(h)
                if allpks.used(h(g))==1
%                     disp('ciao');
                    % get prev vec to compare
                    ref=allpks.pks(:,h(g));
                    % calculate angle b/w vecs
                    sep = min( abs(acosd( dot( curr, ref) / (norm(curr)*norm(ref)))), ...
                        abs(acosd( dot( curr, -ref) / (norm(curr)*norm(-ref)))) );
                    allpks.seps(m,h(g)) = sep;
                    if sep<peaksep
                        sepflag=0;
    %                     break 
                    end
                end
            end            
            % if separated from previous peaks, save
            if sepflag==1
                allpks.used(m) = 1;
                pks.pkvecs = cat(2,pks.pkvecs,curr);
                pks.pkvals = cat(2,pks.pkvals,allpks.vals(m));
                pks.pkfracs = cat(2,pks.pkfracs,allpks.nvals(m));
            end
        end           
    end       
end

allpks.pks = allpks.pks';
pks.pkvecs = pks.pkvecs';

end