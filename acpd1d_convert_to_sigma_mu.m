function [Conductivity,ProductSigmaMu] = acpd1d_convert_to_sigma_mu(Conductivity,Alpha,Eta,Mu0Sigma0,PerturbationParameter)

% note that Mu0Sigma0 is the product of the reference values of mu and
% sigma.



% calculate the Conductivity and permeability then save the result:
val = 1/(1 - PerturbationParameter*Alpha(end)); % value of  1 + epsilon *alpha(0).
Conductivity.RelativeCoef = PerturbationParameter*Alpha*val+1;
Conductivity.CoefValue = Conductivity.RefCoef*Conductivity.RelativeCoef;

ProductSigmaMu = Conductivity; % product of sigma and mu.
ProductSigmaMu.RefCoef = Mu0Sigma0; 

ProductSigmaMu.RelativeCoef = PerturbationParameter*Eta*val+1;
ProductSigmaMu.CoefValue = (ProductSigmaMu.RefCoef)*(ProductSigmaMu.RelativeCoef);
