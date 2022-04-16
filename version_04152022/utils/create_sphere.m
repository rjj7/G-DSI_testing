function [ diff_disp_sphere ]  = create_sphere( sizesphere )

[sx, sy, sz]         = sphere(sizesphere);
diff_disp_sphere     = surf2patch(sx,sy,sz,'triangles');

end