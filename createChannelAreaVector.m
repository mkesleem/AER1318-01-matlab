function [sVec, dsdxVec] = createChannelAreaVector(xVec, xDomain)
%createChannelAreaVector Create vectors for channel area
%   Call functions to fill channel area vector

% Find number of nodes
thisM = length(xVec);
% Create empty matrices
sVec = zeros(thisM,1);
dsdxVec = zeros(thisM,1);

% Calculate values for interior nodes
for idx = 1:thisM
    [thisVal, thisDerivative] = calcChannelArea(xVec(idx));
    sVec(idx) = thisVal;
    dsdxVec(idx) = thisDerivative;
end

end
