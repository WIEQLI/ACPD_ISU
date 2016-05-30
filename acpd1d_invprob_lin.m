function [Conductivity, ProductSigmaMu] = acpd1d_invprob_lin(ParameterFile,DataFileName,sigma0,Mu0)
% solve the linearized inverse problem for the AC potential drop. 
% Input: 
%     ParameterFile: text file containing input parameters
%     DataFileName: name of the potential drop data
%     sigma0, Mu0: reference values for conductivity and relative
%     permeability. These values can be estimated with a homogeneous model.
%     
% Output: 
%    Conductivity and Product of Conductivity and Permeability, two classes.
% =========================================================================

MU0 = 12.5663706143592e-07; % permeabliity of the free space

% % ---- load input parameters: 
[rho11,rho12,rho21,rho22,PerturbationParameter,~,Depth,RegPar] = acpd1d_invprob_lin_get_parameters(ParameterFile);
% PerturbationParameter EPSILON is an estimate of the level of perturbation
% of the coefficient from the reference value for both conductivity and permeability. 

mu0 = Mu0*MU0; % reference value of mu. 

% create the objects for conductivity and permeability: 
Conductivity = PwLinCoefficient(Depth,0*Depth,sigma0); % the values of the coefficients are to be reconstructed

% load the potential drop data file, including frequencies:
dat = dlmread(DataFileName);
Freq = dat(:,1); % set of frequencies
PotentialDrop = dat(:,2) + dat(:,3)*1i; % the second column is the real part, the third column is the imaginary part of the potential drop.


% ==== MAIN PART: solve the linearized optimization problem: call the core function

[Alpha, Eta] = acpd1d_invprob_lin_core(Freq,PotentialDrop,sigma0,mu0,...
              Conductivity, rho11,rho12,rho21,rho22,PerturbationParameter,RegPar);

% ==== END MAIN PART. 

% calculate the conductivity sigma and permeability mu and their product  from Alpha and Eta:
[Conductivity,ProductSigmaMu] = acpd1d_convert_to_sigma_mu(Conductivity,Alpha,Eta,mu0*sigma0,PerturbationParameter);




