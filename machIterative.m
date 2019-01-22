function [machNew] = machIterative(sRatio, mach0)
%machIterative Solve equation 3.45 in the course text iteratively for the
%Mach number

GAMMA = 1.4; % Gamma is 1.4 for air

machFunc = @(m) (2/(GAMMA+1)...
    * (1+(GAMMA-1)/2 * m^2))^( (GAMMA+1)/2/(GAMMA-1) )...
    - m * sRatio;

machDeriv = @(m) m*(2/(GAMMA+1)*...
    (1+(GAMMA-1)/2*m^2))^((3-GAMMA)/2/(GAMMA-1)) - sRatio;

error = 1000;
machOld = mach0;

while error > 5e-14
    machNew = machOld - machFunc(machOld) / machDeriv(machOld);
    error = abs(machNew - machOld) / abs(machOld);
    machOld = machNew;
end

end

