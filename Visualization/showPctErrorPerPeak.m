function showPctErrorPerPeak(IArray,tArray,testToShow,peakNames,...
    paramToTest,paramRange)
%SHOWPCTERRORPERPEAK Summary of this function goes here
%   Detailed explanation goes here

if strcmp(paramToTest,'imageSizeCell')
    legendStr = arrayfun(@(x) ['N = ' num2str(x) ' px'],paramRange,...
        'UniformOutput',false);
else
    legendStr = arrayfun(@(x) [paramToTest ' = ' num2str(x)],paramRange,...
        'UniformOutput',false);
end

intDiff = abs(IArray(:,:,testToShow+1)...
    -IArray(:,:,testToShow));
pctE = 100*intDiff./IArray(:,:,testToShow);
nPlots = numel(peakNames);

figure;
colorList = jet(nPlots).*0.8;
for iPlot = 1:nPlots
    semilogy(tArray,pctE(iPlot,:),...
        'Color',colorList(iPlot,:),'LineWidth',1.5)
    hold on
end
xlabel('Thickness (nm)')
ylabel('(I_{k+1} - I_{k}) / I_{k} (%)')
title(legendStr{testToShow})
legend(peakNames)

end

