clear
close all

% Problem Parameters
p01 = 100e3; % kPa
T01 = 300; % K
GAMMA = 1.4;
R = 287; % Nm/kg/K
criticalArea = 0.8;

% Grid parameters
x0 = 0.0;
x1 = 10.0;
nodeCount = 53;

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
    thisMach = machIterative(thisRatio, thisMach);
    thisP = p01*(1+(GAMMA-1)/2 * thisMach^2)^(-(GAMMA/(GAMMA-1)));
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
rho0 = p01 / R / T01;
qMatExact = zeros(length(xVec),3);
for idx = 1:length(xVec)
    thisMach = machVec(idx);
    thisP = pVec(idx);
    thisRho = rho0 * (1+(GAMMA-1)/2*thisMach^2)^(-1/(GAMMA-1));
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
str1 = "subsonicExactSol";
str2 = sprintf("%d",nodeCount - 2);
str3 = ".mat";
fileNameString = str1 + str2 + str3;
save(fileNameString,"qMatExact", "xVecExact",...
    "machVecExact", "pVecExact");
% Save boundary conditions
qBoundaryMat = zeros(2,3);
qBoundaryMat(1,:) = qMatExact(1,:);
qBoundaryMat(end,:) = qMatExact(end,:);
save("subsonicBoundaryMat.mat","qBoundaryMat");
