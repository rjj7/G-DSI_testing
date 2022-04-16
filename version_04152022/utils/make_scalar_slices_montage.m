function [fig ] = make_scalar_slices_montage( DATAIN, DATAOUT, zslice )

maskslice = squeeze(DATAIN.mask(:,:,zslice));
bkgmask = maskslice==0;

rslice1 = squeeze(DATAIN.lowb(:,:,zslice));
rslice1(bkgmask==1)=0;
rslice2 = squeeze(DATAOUT.scalarmaps.rtop1(:,:,zslice));
rslice2(bkgmask==1)=0;
rslice3 = abs(squeeze(DATAOUT.scalarmaps.std_disp_sph_mean(:,:,zslice)));
rslice3(bkgmask==1)=0;
rslice4 = squeeze(DATAOUT.scalarmaps.kurtosis(:,:,zslice));
rslice4(bkgmask==1)=0;
rslice5 = squeeze(DATAOUT.scalarmaps.isokldiv(:,:,zslice));
rslice5(bkgmask==1)=0;



figH = 5;
figW = 15;
fig = figure('units','inches','position',[0.25 0.25 figW figH],...
    'color','w','InvertHardcopy','off');

[positions] = subplot_pos(9, 4, 0.5, 0.5, 0.5, 0.5, 5, 1, 0.2, 0.2);

ax0 = axes('position',positions{1,1});
imshow(rot90(rslice1),[0 prctile(rslice1,99,'all')]);
title('b=0'); 
colorbar('southoutside');
ax0.Title.FontSize = 16;

ax1 = axes('position',positions{2,1});
imshow(rot90(rslice2),[0 0.25]);
title('RTOP'); 
colorbar('southoutside');
ax1.Title.FontSize = 16;

ax2 = axes('position',positions{3,1});
imshow(rot90(rslice3),[0 0.03]);
title('Mean Radial Std. Dev.'); 
colorbar('southoutside');
ax2.Title.FontSize = 16; 

ax3 = axes('position',positions{4,1});
imshow(rot90(rslice4),[2.2 2.85]);
title('Mean kurtosis');
colorbar('southoutside');
ax3.Title.FontSize = 16; 

ax4 = axes('position',positions{5,1});
imshow(rot90(rslice5),[0.2 1]);
title('AnIso KLDiv');
colorbar('southoutside');
ax4.Title.FontSize = 16;


end
