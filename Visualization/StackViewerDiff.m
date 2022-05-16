function [] = StackViewerDiff(imageStack,tArray)
%STACKVIEWER_V2 User can browse images in a M X N X P stack along the
%P dimension (assumed to be the stack dimension). A slider bar is given at
%the bottom for scrolling through images.
%   
% To fix:
% - Allow window size to be specified (could even just be small, medium,
% large... ideally would choose based on computer screen size)
% - Meaningful x and y axes
% - Colorbar/scale

% Building the main figure window 
fig = figure('NumberTitle','off','Name','Image Stack Viewer');
set(fig,'Position',[0 0 500 600]);
set(fig,'Color',[0.9 0.9 0.9]);
set(fig,'Units','pixels');

% homemade percentile function
pctl = @(v,p) interp1(linspace(0.5/length(v), 1-0.5/length(v), length(v))',...
    sort(v), p*0.01, 'spline');
clims = [0 pctl(imageStack(:),99)];

% Setup the image display axis
im = imagesc(imageStack(:,:,1));
axis off
set(gca,'Units','pixels');
set(gca,'Position',[50 100 350 350]);
% caxis([0 5000]);
caxis(clims)
colormap(parula(1024))
axis ('square')
title(['I(q_x,q_y) at t = ' num2str(tArray(1),3) ' nm'])
% uicontrol('Parent',fig,'Style','text','String','1',...
%           'Position',[250 20 20 20]);

% Create the slider bar
dim3 = size(imageStack,3);
frameslider = uicontrol('Parent',fig,'Style','slider',...
               'Min',1,'Max',dim3,'Value',1,...
               'SliderStep',[1/dim3 1/dim3],...
               'Position',[50 50 350 20],...
               'BackgroundColor',[1 1 1],...
               'Callback',@frameCallback);
% set(frameslider,'BackgroundColor',[1 1 1])
% set(frameslider,'Callback',@frameCallback);

% Callback function for interactive update of image
function[] = frameCallback(source,~)
    pos = get(source,'value');
    pos = int16(pos);
    set(im,'CData',imageStack(:,:,pos));
    title(['I(q_x,q_y) at t = ' num2str(tArray(pos),3) ' nm'])
end

end

