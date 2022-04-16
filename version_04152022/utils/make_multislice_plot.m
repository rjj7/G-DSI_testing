function [fig] = make_multislice_plot( data, slicedim, slice, frames, titles, showit)
% 
% [fig] = make_multislice_plot( data, slicedim, slice, frames, titles, showit)
% 
% Make figure plotting one slice for different data frames (4th dim)
% 
%   - titles is optional (can be cell array with fig titles)
%   - showit is optional (true will display figure + is def.)
% 
%   - fig returns figure handle
%   

if nargin<5
    titles = [];
end
if nargin<6
    showit = 1;
end


fig = figure('color','w','InvertHardcopy','off','visible','off',...
    'position',[51 518 1868 433]);

if slicedim==1
    data = squeeze(data(slice,:,:,:));
elseif slicedim==2
    data = squeeze(data(:,slice,:,:));
elseif slicedim==3
    data = squeeze(data(:,:,slice,:));
end

nframes = length(frames);
for id=1:nframes
    frame = frames(id);
    tmp = data(:,:,frame);
    caxislims = [0 prctile(tmp(tmp>0),95)];
    
    subplot(1,nframes,id);
    imshow(tmp,caxislims);
    if ~isempty(titles)
        title(titles{id});
    end

end

if showit
    fig.Visible='on';
end

end
