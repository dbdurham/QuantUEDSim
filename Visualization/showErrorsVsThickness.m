function showErrorsVsThickness(tArray,R,MPctE,paramToTest,paramRange)
%SHOWERRORSVSTHICKNESS Summary of this function goes here
%   Detailed explanation goes here

if strcmp(paramToTest,'imageSizeCell')
    legendStr = arrayfun(@(x) ['N = ' num2str(2*x) ' px'],paramRange,...
        'UniformOutput',false);
else
    legendStr = arrayfun(@(x) [paramToTest ' = ' num2str(x)],paramRange,...
        'UniformOutput',false);
end

figure('Position',[100,100,900,300]);
subplot(1,2,1)
semilogy(tArray,R,'LineWidth',1.5)
xlabel('Thickness (nm)')
ylabel('R (%)')
legend(legendStr{1:end-1})
subplot(1,2,2)
semilogy(tArray,MPctE,'LineWidth',1.5)
xlabel('Thickness (nm)')
ylabel('Mean % Error')
legend(legendStr{1:end-1})

end

