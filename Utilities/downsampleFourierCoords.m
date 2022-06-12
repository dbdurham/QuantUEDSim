function [qxaStore,qyaStore,storeMask] = downsampleFourierCoords(...
    qxa,qya,downSampFac,varargin)
%DOWNSAMPLEFOURIERCOORDS Summary of this function goes here
%   Detailed explanation goes here

if nargin>3
    storeMask = varargin{1};
else
    storeMask = ones(size(qxa));
end

dqx = qxa(2,1)-qxa(1,1);
dqy = qya(1,2)-qya(1,1);
pxa = round(qxa/dqx); % pixels away from center
pya = round(qya/dqy);
if downSampFac > 1
    storeMask = storeMask & ...
         mod(pxa,downSampFac) == 0 & ...
         mod(pya,downSampFac) == 0;
end
storeLenX = sum(storeMask(:,1));
storeLenY = sum(storeMask(1,:));
qxaStore = reshape(qxa(storeMask),...
    [storeLenY,storeLenX]);
qyaStore = reshape(qya(storeMask),...
    [storeLenY,storeLenX]);


end

