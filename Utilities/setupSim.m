function sDiff = setupSim(simType,options)
%setupSim Wrapper for setup functions
%   
if nargin == 1
    options = struct;
end

switch simType
    case 'Kinematical'
        sDiff = setupSimKin(options);
    case 'Bloch Waves'
        sDiff = setupSimBW(options);
    case 'Multislice'
        sDiff = setupSimMS(options);
end

end

