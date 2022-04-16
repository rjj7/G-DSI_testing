function [ff]=plot_ralpha_contour(ralpha, diff_disp_sphere)    

% For a voxel, make pdf contour glyphs
    pdf_actor = diff_disp_sphere;
    pdf_actor.vertices = diff_disp_sphere.vertices .* repmat(ralpha, [1, 3]); % scale radial distance
    pdf_actor.facevertexcdata = ralpha; % change color data to represent pdf values
    if min(pdf_actor.faces(:))==0, pdf_actor.faces=pdf_actor.faces+1; end

    % display pdf contour
    ff=figure('color','w','position',[288 415 855 460]);
    ff.InvertHardcopy='off';
%                 subplot(1,2,1);
    h = patch(pdf_actor);
    view(180,0); %view(90, 90);
    lighting gouraud; shading faceted; camlight
    set(h, 'EdgeColor', 'none');
    % colormap hsv; %colorbar;
    caxis([min(ralpha(:)), max(ralpha(:))]);
    axis equal, axis off, axis tight
%                 title(['pdf contour ' sprintf('%.03f',r0) 'xMDD']);
  
