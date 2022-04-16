function [ diff_disp_sphere ]  = create_sphere_nneg( sizesphere )

[sx, sy, sz]         = sphere(sizesphere);

sx = sx(((sizesphere+1)/2)+1:end,:);
sy = sy(((sizesphere+1)/2)+1:end,:);
sz = sz(((sizesphere+1)/2)+1:end,:);

diff_disp_sphere     = surf2patch(sx,sy,sz);

end