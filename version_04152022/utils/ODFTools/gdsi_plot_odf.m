function [ f ] = gdsi_plot_odf( fvc, ODF, ODFpeaks )
% [ f ] = gdsi_plot_odf( fvc, ODFpoints[, ODFpeaks] )
% -- expects ODFpeaks to be vecs size [nx3]

if nargin<3
    ODFpeaks=[];
end

% For a voxel, make pdf contour glyphs
pdf_actor = fvc;
if min(pdf_actor.faces(:))==0
    pdf_actor.faces = pdf_actor.faces+1;
end
ODFpoints = fvc.vertices .* repmat(ODF', [1, 3]);
pdf_actor.vertices = ODFpoints; % scale radial distance
s = vecnorm(ODFpoints');
ss = ODFpoints./repmat(s',1,3);
pdf_actor.facevertexcdata = abs(ss); % change color data to represent pdf values

% display pdf contour
f=figure('color','w','position',[288 415 500 500]);
f.InvertHardcopy='off';

h = patch(pdf_actor);
view(180,0);
lighting gouraud; 
shading faceted; 
camlight
set(h, 'EdgeColor', 'none');
% colormap hsv; %colorbar;
% caxis([min(ODF(:)), max(ODF(:))]);
axis equal, axis off, axis tight
%                 title(['pdf contour ' sprintf('%.03f',r0) 'xMDD']);

if ~isempty(ODFpeaks)
    hold on;
    pkcolors = {'k','c','m'};
    nmax=min(3,size(ODFpeaks,1));
    for nn=1:nmax
        ccolor=pkcolors{nn};
        currvec = ODFpeaks(nn,:);
        px = 1.25*[currvec(1) -currvec(1)];
        py = 1.25*[currvec(2) -currvec(2)];
        pz = 1.25*[currvec(3) -currvec(3)];
        plot3(px,py,pz,'-','LineWidth',3,'Color',ccolor);
    end
    
end