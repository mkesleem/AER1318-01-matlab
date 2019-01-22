function [pNew] = pressureIterative(p0, paramStruct)
%pressureIterative Solve iteratively for p2 for the shock-tube problem

% Extract useful values
pR = paramStruct.pR;
pL = paramStruct.pL;
aR = paramStruct.aR;
aL = paramStruct.aL;
gamma = paramStruct.gamma;
alpha = paramStruct.alpha;

% Define functions for Newton's method
pFunc = @(p) sqrt(2/gamma/(gamma-1)) * (p-1)/sqrt(1+alpha*p)...
    -2/(gamma-1) * aL/aR * (1 - (pR/pL * p)^((gamma-1)/2/gamma));
pDeriv = @(p) 1/2*sqrt(2/gamma/(gamma-1))*((3+alpha*p)/(1+alpha*p)^(3/2))...
    + 1/gamma * aL/aR * pR/pL * (pR/pL*p)^(-(gamma+1)/2/gamma);

% Implement Newton's method
error = 1000;
pOld = p0;

while error > 5e-14
    pNew = pOld - pFunc(pOld) / pDeriv(pOld);
    error = abs(pNew - pOld) / abs(pOld);
    pOld = pNew;
end

end

