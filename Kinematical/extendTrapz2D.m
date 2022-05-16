function [Itrapz,Isimp,xList,yList,wList] = extendTrapz2D(func, x1,x2,y1,y2,I,n)
%EXTENDTRAPZ Apply extended trapezoidal rule to iteratively integrate a
%function using a square sampling grid
%   func -- function to integrate, set up as I = func(x,y)
%   x1, x2 -- limits in x (column direction)
%   y1, y2 -- limits in y (row direction)
%   I -- prior integral value
%   n -- current iteration; n = 1 is the first corner step
% 
%   Returns:
%   Itrapz -- new integral value via trapezoidal rule (feed into next
%   iteration)
%   Isimp -- new integral value via Simpson rule (calculated from new and
%   previous trapezoidal rule iterations)
% 

if n == 1
    dely = (y2-y1);
    delx = (x2-x1);
    Itrapz = 1/4.*(func(x1,y1)+func(x1,y2)+func(x2,y1)+func(x2,y2))*delx*dely; % sum of corners
    Isimp = 0; % Insufficient sampling to apply Simpson's rule yet
    xList = [x1,x1,x2,x2];
    yList = [y1,y2,y1,y2];
    wList = [1/4,1/4,1/4,1/4];
else
    % generate list of x,y coordinates at which to sample, plus
    % corresponding weights
    
    % define useful variables
    Nrows = 2 + sum(2.^((2:n)-2)); % number of rows in sampling grid
    Ncols = Nrows; % number of columns in sampling grid
    Nsamps = Nrows^2; % number of total samples in sampling grid
    NrowsPrev = 2 + sum(2.^((2:n-1)-2));
    NsampsPrev = NrowsPrev^2; % number of samples already computed in previous grid
    NsampsNew = Nsamps-NsampsPrev; % number of new samples to compute
    NcolsNew = Nrows-NrowsPrev; % number of new columns
    % sampling grid spacings
    dely = (y2-y1)/(Nrows-1);
    delx = (x2-x1)/(Ncols-1);
    
    % For reference: position and weights of samples
    xList = zeros(NsampsNew,1);
    yList = zeros(NsampsNew,1);
    wList = zeros(NsampsNew,1);
    
    Inew = zeros(size(I));
    iSampTot = 1;
    for iRow = 1:Nrows
        y = y1+(iRow-1)*dely;
        if mod(iRow,2) == 0 % even row
            for iSamp = 1:Ncols % a new sample in each column
                x = x1+(iSamp-1)*delx;
                if (iSamp == 1 || iSamp == Ncols)
                    Inew = Inew + 0.5*func(x,y)*delx*dely;
                    wList(iSampTot) = 0.5;
                else
                    Inew = Inew + func(x,y)*delx*dely;
                    wList(iSampTot) = 1;
                end
                xList(iSampTot) = x;
                yList(iSampTot) = y;
                iSampTot = iSampTot + 1;
            end
        else % odd row
            for iSamp = 1:NcolsNew % new samples only in new columns
                x = x1 + delx + (iSamp-1)*2*delx;
                if (iRow == 1 || iRow == Nrows)
                    Inew = Inew + 0.5*func(x,y)*delx*dely;
                    wList(iSampTot) = 0.5;
                else
                    Inew = Inew + func(x,y)*delx*dely;
                    wList(iSampTot) = 1;
                end
                xList(iSampTot) = x;
                yList(iSampTot) = y;
                iSampTot = iSampTot + 1;
            end
        end
    end
    
    Itrapz = Inew + (1/4)*I; % Add old contribution to the integral
    Isimp = (16/15)*Itrapz - (1/15)*I; % Simpson rule integral
end

end


