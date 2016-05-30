function [Alpha, Eta,fval] = acpd1d_invprob_lin_core(Freq,PotentialDrop,sigma0,mu0,Conductivity, rho11,rho12,rho21,rho22,PerturbationParameter,RegPar)
% this is the core function for solving linearized ACPD problem. Both
% recursive and non-recursive algorithms use this function. 
% Input: Freq: a vector of frequencies
% PotentialDrop: a vector of potential drop (complex numbers) with the same number of elements as in Freq. 
% sigma0, mu0: reference values of conductivity and permeability, respectively.
% Conductivity: an object of class "Coefficient" (or its subclasses), just
% use the grid to calcuate the coefficient matrices
% rho11,rho12,rho21,rho22: distances between the pins
% PerturbationParameter: (epsilon): estimated level of perturbation of the coefficients from the reference values
% OUTPUT: 
%    Alpha(y) = alpha(y)/[1 + epsilon*alpha(0)]; 
%    Eta(y) = eta(y)/[1 + epsilon*alpha(0)]
% 
% Note: from Alpha and Eta, we will convert to alpha and eta. 
% Created March 25, 2016.


Nz = length(Conductivity.Depth); % number of grid points in depth


% calculate the zero order term and the right hand side vector for the linearized model: 
D0 = acpd1d_function_D0(Freq,sigma0,mu0,rho11,rho12,rho21,rho22); % the zero order term
RHS = (PotentialDrop - D0(:))/PerturbationParameter; % the right hand side vector of the linearized problem

% Compute the coefficient matrix of the linearized inverse problem: 
[A,B] = Conductivity.CoefficientMatrix(Freq,sigma0,mu0,rho11,rho12,rho21,rho22);   

% %+++ ---test with the exact solution (for simulated data, if the exact solution is known): 
%     Alpha = (ConductivityInit.RelativeCoef-1)/(ConductivityInit.RelativeCoef(Nz));
%     Eta = (PermeabilityInit.RelativeCoef-1)/(ConductivityInit.RelativeCoef(Nz));
% %     Alpha = (ConductivityInit.RelativeCoef-1);
% %     Eta = (PermeabilityInit.RelativeCoef-1);
%     PD = D0 + A*Alpha + B*Eta;        
%     max(abs(PD - PDData))
% % +++ ---end test!                
        
% scale the coefficient matrix and the right hand side:
RHS = RHS*sigma0; RHS = RHS(:);
A = A*sigma0; B = B*sigma0;


%--solve the linear constrained opt. problem as the quadratic programming
CoefMat = [[real(A), real(B)]; [imag(A), imag(B)]]; % combine both real and imaginary parts
RHS = [real(RHS); imag(RHS)];

RHS = (CoefMat'*RHS);
CoefMatrix = (RegPar*eye(2*Nz) + CoefMat'*CoefMat);

% Sol = CoefMatrix\RHS


% Optimization approach:  
%constraints: 
ConstraintRHS = zeros(4*Nz-6,1);
ZeroMat1 = zeros(Nz-1,Nz); 
ConstMat1 = ZeroMat1;

for n = 1:Nz-1
    ConstMat1(n,n) = -1;
    ConstMat1(n,n+1) = 1;
end
ConstMatrix1 = [[ConstMat1, ZeroMat1]; [ZeroMat1, ConstMat1]];

ZeroMat2 = zeros(Nz-2,Nz);
ConstMat2 = ZeroMat2;

for n = 1:Nz-2
    ConstMat2(n,n) = 1;
    ConstMat2(n,n+1) = -2;
    ConstMat2(n,n+2) = 1;
end
ConstMatrix2 = [[ConstMat2, ZeroMat2]; [ZeroMat2, ConstMat2]];
ConstraintMatrix = [ConstMatrix1; ConstMatrix2];


% ----solving the optimization problem:
% options = optimoptions('quadprog','display','iter','algorithm','interior-point-convex','tolFun',1e-15);
options = optimset('display','iter','algorithm','interior-point-convex','tolFun',1e-15);
Aeq = []; beq = []; % no equality constraints
LowerBound = -3*ones(2*Nz,1); 
UpperBound = 3*ones(2*Nz,1);

[Sol,fval] = quadprog(CoefMatrix,-RHS,ConstraintMatrix,ConstraintRHS,Aeq,beq,LowerBound,UpperBound,ones(2*Nz,1)/2,options); % no equality/bound constraints. 

plot(fval);
% Sol = cgm_uncon_linear(CoefMatrix,RHS,zeros(2*Nz,1),1e-14);
cond(CoefMatrix)

% % end of optimization approach

% -----extract the estimated coefficients: 
Alpha = Sol(1:Nz); 
Eta = Sol(Nz+1:2*Nz);
