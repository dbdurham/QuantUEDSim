function sDiff = setupSim(simType,varargin)
%setupSim Wrapper for setup functions
%   
if nargin > 1
    options = varargin{1};
else
    options = struct;
end

switch simType
    case 'Kinematical'
        sDiff = setupSimKin(varargin);
    case 'Bloch Waves'
        sDiff = setupSimBW(options);
    case 'Multislice'
        sDiff = setupSimMS(options);
end

end

