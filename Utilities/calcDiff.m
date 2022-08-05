function IDiff = calcDiff(theta1,theta2,numUC,sDiff)
%CALCDIFF Wrapper function for all diffraction pattern computations
%   

switch sDiff.simType
    case 'Kinematical'
        IDiff = calcDiffKin(theta1,theta2,numUC,sDiff);
    case 'Bloch Waves'
        IDiff = calcDiffBW(theta1,theta2,numUC,sDiff);
    case 'Multislice'
        IDiff = calcDiffMS(theta1,theta2,numUC,sDiff);
end

end

