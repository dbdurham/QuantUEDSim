function [y,Y] = multiExpConvGauss(c,x,n)
% MULTIEXPCONVGAUSS generates mutiple exponential functions convolved with 
% an identical gaussian kernel. Used for modeling time traces with gaussian
% instrument response.
    % n is the number of functions
    % c is the parameters: gaussian sigma, time zero, exp amplitudes (n),
    % exp time constants (n), exponential offsets (n)
    % x is the x position data
    % y is the output as a flattened vector
    % Y is the output as an n by length(x) matrix
    
    
    % use 1000 points covering twice the range to compute conv(gauss,exp)
    xConv = linspace(x(1)-0.5*(x(end)-x(1)),x(end)+0.5*(x(end)-x(1)),1000);
    dx = 2*(x(end)-x(1))/1000;
    XConv = repmat(xConv,[n 1]);
    xKernel = -4*c(1):dx:4*c(1);
    yKernel = (1/(c(1)*sqrt(2*pi)))*exp(-0.5*(xKernel./(c(1))).^2);
    X = repmat(x,[n 1]);
    CAmp = repmat(c(3:2+n)',[1 length(xConv)]);
    CTau = repmat(c(3+n:2+2*n)',[1 length(xConv)]);
    COffset = repmat(c(3+2*n:end)',[1 length(xConv)]);
    
    YExp = (XConv>c(2)).*CAmp ...
    .*(exp(-(XConv-c(2))./CTau) - 1)...
    + COffset;
% %   DBD: attempt to avoid for loops
%     YKernel = [zeros(ceil((n-1)/2),length(yKernel));...
%         yKernel;...
%         zeros(floor((n-1)/2),length(yKernel))];
%     YConv = conv2(YExp,YKernel,'same');

    Y = zeros(size(X));
    for ii = 1:n
        yConv = conv(YExp(ii,:),yKernel,'same').*dx;
        Y(ii,:) = interp1(xConv,yConv,x);
    end
        
    y = Y(:);
end


