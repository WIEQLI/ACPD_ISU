clear all
N = 6;
depth = linspace(-1,0,N);
coefvalue = rand(1,N);
refcoef = 2;

permittivity = Coefficient(depth,coefvalue,refcoef);
depth2 = permittivity.Depth; 
CoefValue = permittivity.CoefValue;


permittivity.plot();


% subclass: piecewise linear coefficient: 
perm2 = PwLinCoefficient(depth,coefvalue,refcoef);
perm2.plot2();


% test the coefficient matrix calculation: 
Freq = linspace(10,100,5); sigma0 = 1; mu0 =  1; rho11 = 0.0015; rho12 = 0.003; rho21 = rho12; rho22 = rho11;
[A,B] = CoefficientMatrix(perm2,Freq,sigma0,mu0,rho11,rho12,rho21,rho22);


% interpolate data: 
z = linspace(-1,0,2*N-1);
perm2_2 = perm2.Interpolation(z);
perm2_2.plot();

% subclass: piecewise constant coefficient:
perm3 = PwConstCoefficient(depth,coefvalue,refcoef);
perm3.plot();