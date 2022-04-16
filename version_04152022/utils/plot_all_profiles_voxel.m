function [ f ] = plot_all_profiles_voxel( x, data, vertices )

f = figure; hold on;
f.Color='w';
f.InvertHardcopy='off';
f.Position=[1228 96 619 557];

nverts = size(data,2);
if nverts~=length(vertices)
    vertices = vertices(:,1:nverts);
end
vertices = abs(vertices);

for vv=1:nverts
    tmp = data(:,vv);
    plot(x,tmp,'Color',vertices(:,vv));
end

xlabel('Displacement (microns)'); 
ylabel('Probability density (a.u.)');


