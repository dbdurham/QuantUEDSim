function [m,sm] = linRegMLENoInt(x,y,s)
%LINREGMLENOINT Computes the best fit slope and error for a 
% no-intercept linear regression using maximum likelihood estimation
% likelihood estimation 
%   x - independent variable values
%   y - dependent variable values
%   s - Gaussian standard error values

% y = mx 
% chi^2 = sum( ((y - m*x)./s).^2)
% Minimize by setting derivative w.r.s. m to zero

% Summation shorthands
Sxy = sum((x.*y) ./ s.^2);
Sxx = sum(x.^2 ./ s.^2);

m = Sxy/Sxx;
sm = sqrt(1/Sxx);

end

