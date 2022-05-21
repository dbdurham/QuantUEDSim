function showIvtVsParam(IArray,tArray,peakNames,paramToTest,paramRange)
%SHOWIVTVSPARAM Summary of this function goes here
%   Detailed explanation goes here

if strcmp(paramToTest,'imageSizeCell')
    legendStr = arrayfun(@(x) ['N = ' num2str(x) ' px'],paramRange,...
        'UniformOutput',false);
else
    legendStr = arrayfun(@(x) [paramToTest ' = ' num2str(x)],paramRange,...
        'UniformOutput',false);
end
figure
nPlots = numel(peakNames);
for iPlot = 1:nPlots
    subplot(2,ceil(nPlots/2),iPlot)
    plot(tArray,squeeze(IArray(iPlot,:,:)),...
        '-o')
    xlabel('Thickness (nm)')
    ylabel('Fraction of diffracted intensity')
    title(peakNames{iPlot})
    if iPlot == nPlots
       lg= legend(legendStr);
    end
end

% 
% set(gcf,'color','white','position',[50 50 850 400]);
% set(lg,'position',[0.85 0.35 0.1 0.1]);

end

