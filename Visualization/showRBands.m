function showRBands(RBands,tBands,paramRange,paramToTest)
%SHOWRBANDS Summary of this function goes here
%   Detailed explanation goes here
nBands = size(RBands,1);
colorList = 0.8*jet(nBands);
figure;
for iBand = 1:nBands
    if nBands == 1
        plot(paramRange,RBands,'k.-',...
        'MarkerSize',16,'LineWidth',1.5)
    else
        semilogy(paramRange,RBands(iBand,:),'.-',...
            'MarkerSize',16,'LineWidth',1.5,...
            'Color',colorList(iBand,:))
    end
    hold on
end
plot([0 paramRange(end)*1.1],[1 1],'k--','LineWidth',1.5)
if strcmp(paramToTest,'imageSizeCell')
    xlabel('Image width (px)')
else
    xlabel(paramToTest)
end
ylabel('R (%)')

legendStr = cell(nBands,1);
for iBand = 1:nBands
    legendStr{iBand} = ...
        ['t = ' num2str(tBands(iBand,1),2) ...
        ' - ' num2str(tBands(iBand,2),2) ' nm'];
end
legend(legendStr)

end

