clear
close all

% Problem Parameters
p01 = 100e3; % kPa
T01 = 300; % K
GAMMA = 1.4;
R = 287; % Nm/kg/K
criticalArea = 1;

% Grid parameters
x0 = 0.0;
x1 = 10.0;
nodeCount = 201;

% Create vector x
xVec = linspace(x0, x1, nodeCount);
xDomain = [x0, x1];
deltaX = (x1 - x0) / (nodeCount-1);

% Calculate channel area
[sVec, dsdxVec] =...
    createChannelAreaVector(xVec, xDomain);

pVec = zeros(size(xVec));
machVec = zeros(size(xVec));
thisMach = 0.1;

for idx = 1:length(xVec)
    thisRatio = sVec(idx) / criticalArea;
    thisMach = machIterative(thisRatio, thisMach+0.001);
    thisP = p01*(1+(GAMMA-1)/2 * thisMach^2)^(-(GAMMA/(GAMMA-1)));
    
    if abs((xVec(idx) - 7)) < 5e-15
        machL = thisMach;
        pL = thisP;
        p0L = p01;
        T0L = T01;
        criticalAreaL = criticalArea;
        % Mach number
        machRtop = 2 + (GAMMA-1) * machL^2;
        machRbot = 2*GAMMA*machL^2 - (GAMMA-1);
        machR = sqrt(machRtop / machRbot);
        % Pressure
        pRtop = 2 * GAMMA * machL^2 - (GAMMA-1);
        pRbot = GAMMA + 1;
        pR = pRtop / pRbot * pL;
        % Ambient pressure
        p0Rtop = (GAMMA+1)/2 * machL^2 / (1 + (GAMMA-1)/2*machL^2);
        p0Rtop = p0Rtop^(GAMMA/(GAMMA-1));
        p0Rbot = 2*GAMMA / (GAMMA+1) * machL^2 - (GAMMA-1)/(GAMMA+1);
        p0Rbot = p0Rbot^(1/(GAMMA-1));
        p0R = p0Rtop / p0Rbot * p0L;
        % Ambient Temperature
        T0R = T0L;
        
        % Ideal gas relations
        rho0L = p0L / R / T0L;
        rho0R = p0R / R / T0R;
        a0L = sqrt(GAMMA*p0L/rho0L);
        a0R = sqrt(GAMMA*p0R/rho0R);
        aRhoL = rho0L * a0L * (2/(GAMMA+1))^((GAMMA+1)/2/(GAMMA-1));
        aRhoR = rho0R * a0R * (2/(GAMMA+1))^((GAMMA+1)/2/(GAMMA-1));
        
        % Calculate critical area
        criticalAreaR = criticalAreaL * aRhoL / aRhoR;
        
        % Change value of p0 and S
        p01 = p0R;
        criticalArea = criticalAreaR;
        thisMach = machR;
        thisP = pR;
    end
    
    machVec(idx) = thisMach;
    pVec(idx) = thisP;
end

figure(1)
plot(xVec,pVec);
xlabel("x");
ylabel("pressure (Pa)");

figure(2)
plot(xVec,machVec);
xlabel("x");
ylabel("Mach number");

% Solve for q
qMatExact = zeros(length(xVec),3);
for idx = 1:length(xVec)
    thisMach = machVec(idx);
    thisP = pVec(idx);
    if xVec(idx) <= 7
        thisRho = rho0L * (1+(GAMMA-1)/2*thisMach^2)^(-1/(GAMMA-1));
    else
        thisRho = rho0R * (1+(GAMMA-1)/2*thisMach^2)^(-1/(GAMMA-1));
    end
    thisSound = sqrt(GAMMA*thisP/thisRho);
    thisU = thisMach * thisSound;
    thisE = thisP/(GAMMA-1) + thisRho*thisU^2/2;
    q1 = thisRho;
    q2 = thisRho*thisU;
    q3 = thisE;
    qVec = [q1, q2, q3] * sVec(idx);
    qMatExact(idx,:) = qVec;
end

% Save exact solution
pVecExact = pVec;
xVecExact = xVec;
machVecExact = machVec;
save("transonicExactSol","qMatExact", "xVecExact",...
    "machVecExact", "pVecExact");
% Save boundary conditions
qBoundaryMat = zeros(2,3);
qBoundaryMat(1,:) = qMatExact(1,:);
qBoundaryMat(end,:) = qMatExact(end,:);
save("transonicBoundaryMat.mat","qBoundaryMat");
