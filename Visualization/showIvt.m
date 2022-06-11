function showIvt(IArray,I0Array,tArray,peakNames)
%SHOWIVT Plot primary and diffracted beam intensities vs thickness
%   Detailed explanation goes here

nPeaks = size(IArray,1);

figure;
subplot(2,1,1)
plot(tArray,I0Array.^0.5,'k-','LineWidth',1.5)
xticks([])
ylabel('I_{000}^{1/2}')
xlim([0 tArray(end)])
ylim([0 1])

subplot(2,1,2)
nPlots = nPeaks;
colorList = jet(nPlots).*0.8;
for iPlot = 1:nPlots
    plot(tArray,squeeze(IArray(iPlot,:).^0.5),...
            '-','Color',colorList(iPlot,:),'LineWidth',1.5)
    hold on
end
xlabel('Thickness (nm)')
ylabel('I_{Diff}^{1/2}')
legend(peakNames)
xlim([0 tArray(end)])

end

