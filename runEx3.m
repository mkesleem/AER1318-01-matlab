clear
close all

% Problem Parameters
pL = 10^5;
rhoL = 1;
pR = 10^4;
rhoR = 0.125;
GAMMA = 1.4;
t = 6.1e-3;
alpha = (GAMMA+1)/(GAMMA-1);

% Solve for speed of sound
aL = sqrt(GAMMA*pL/rhoL);
aR = sqrt(GAMMA*pR/rhoR);

% Create a parameter structure
paramStruct.pL = pL;
paramStruct.pR = pR;
paramStruct.aL = aL;
paramStruct.aR = aR;
paramStruct.alpha = alpha;
paramStruct.gamma = GAMMA;

% Solve for p2
p = pressureIterative(1, paramStruct);
p2 = p*pR;

% Solve for rho3 and p3
p3 = p2; % Pressure across contact surface is continous
rho3 = rhoL * (p3 / pL)^(1/GAMMA);

% Solve for speeds
v = 2/(GAMMA-1) * aL * (1-(p3 / pL)^((GAMMA-1)/2/GAMMA));
%fluid velocity on either side of the contact surface must be equal to V
u2 = v;
u3 = v;
c = (p-1)*aR^2/GAMMA/u2;

% Create vector x
xDomain = [0.0, 10.0];
nodeCount = 10000;
xVec = linspace(xDomain(1), xDomain(2), nodeCount);
x0 = (xDomain(1) + xDomain(end))/2;

% Solve for original left quiescent state (M=0, rho=1)
xL = x0 - aL * t;

% Solve for state 5
x5L = xL;
x5R = x0 + (v*(GAMMA+1)/2 - aL)*t;
u5Func = @(x) 2/(GAMMA+1)*((x-x0)/t + aL);
a5Func = @(x) u5Func(x) - (x-x0)/t;
p5Func = @(x) pL * (a5Func(x)/aL)^(2*GAMMA / (GAMMA-1));
rho5Func = @(x) GAMMA * p5Func(x) / (a5Func(x))^2;

% Solve for state 3
x3L = x5R;
x3R = x0 + v*t;
a3 = sqrt(GAMMA * p3 / rho3);
mach3 = u3 / a3;

% Solve for state 2
x2L = x3R;
x2R = x0 + c*t;
rho2 = rhoR * (1+alpha*p)/(alpha + p);
a2 = sqrt(GAMMA*p2/rho2);
mach2 = u2 / a2;

% Plot
machVec = zeros(size(xVec));
rhoVec = zeros(size(xVec));
for idx = 1:length(xVec)
    thisX = xVec(idx);
    
    if thisX < xL
        thisMach = 0;
        thisRho = 1;
    elseif thisX >= x5L && thisX < x5R
        thisU = u5Func(thisX);
        thisA = a5Func(thisX);
        thisRho = rho5Func(thisX);
        thisMach = thisU / thisA;
    elseif thisX >= x3L && thisX < x3R
        thisMach = mach3;
        thisRho = rho3;
    elseif thisX >= x2L && thisX < x2R
        thisMach = mach2;
        thisRho = rho2;
    elseif thisX > x2R
        thisMach = 0;
        thisRho = rhoR;
    end
    
    % store values
    machVec(idx) = thisMach;
    rhoVec(idx) = thisRho;
end

figure(1)
plot(xVec,rhoVec);
xlabel("x");
ylabel("\rho (kg/m^3)");
axis([xDomain(1) xDomain(end) 0 1.05]);

figure(2)
plot(xVec,machVec);
xlabel("x");
ylabel("Mach number");