function [Conductivity, ProductSigmaMu] = acpd1d_invprob_lin_recursive(LinParFile,HomParFile,DataFileName)
% solve the linearized inverse problem for the AC potential drop using a recursive algorithm. 
% Input: 
%     LinParFile: text file containing input parameters for the linearized
%     inverse problem
%     HomParFile: parameter file for the homogeneous inverse problem
%     DataFileName: name of the potential drop data
% Output: 
%    Conductivity and Permeability*Conductivity, two classes.
% =========================================================================

MU0 = 12.5663706143592e-07; % permeabliity of the free space

% % ---- load input parameters: 
[rho11,rho12,rho21,rho22,PerturbationParameter,Nz,Depth,RegPar] = acpd1d_invprob_lin_get_parameters(LinParFile);

% load the potential drop data file, including frequencies:
[Freq,PotentialDrop] = acpd1d_invprob_load_data(DataFileName);
NumFreq = length(Freq); % number of frequencies 


    [Sigma0,Mu0] = acpd1d_invprob_hom(HomParFile,Freq,PotentialDrop);
    Sigma0
    Mu0


% ======= Recursive procedure for estimating the coefficients based on different frequency ranges: 
% step 1: use two largest frequencies, estimate constant coefficients:
NoFreq = 2; % number of frequencies in step 1
freq = Freq(NumFreq - NoFreq+1:NumFreq);
potdrop = PotentialDrop(NumFreq - NoFreq+1:NumFreq);
[sigma0,Mu0] = acpd1d_invprob_hom(HomParFile,freq,potdrop,Sigma0,Mu0); % estimate of the coefficient values near the surface.
sigma0
Mu0  % Mu is the relative permeability

mu0 = Mu0*MU0; % reference value of mu 

% Step 2 -- Iteration: increase the number of frequencies gradually
for iter = 1:3

    depth = Depth(); 
    % create the objects for conductivity and permeability: 
    Conductivity = PwLinCoefficient(depth,0*depth,sigma0); % the values of the coefficients to be reconstructed

    freq = Freq(NumFreq - NoFreq+1:NumFreq);
    potdrop = PotentialDrop(NumFreq - NoFreq+1:NumFreq);

    % ==== MAIN PART: solve the linearized optimization problem: call the core function

    [Alpha, Eta] = acpd1d_invprob_lin_core(freq,potdrop,sigma0,mu0,...
                  Conductivity, rho11,rho12,rho21,rho22,PerturbationParameter,RegPar);

    % ==== END MAIN PART. 
end


% -----extract the estimated coefficients: 
Alpha = Sol(1:Nz); 
Eta = Sol(Nz+1:2*Nz);
[Conductivity,ProductSigmaMu] = acpd1d_convert_to_sigma_mu(Conductivity,Alpha,Eta,mu0*sigma0,PerturbationParameter);

