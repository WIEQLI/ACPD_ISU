function run_acpd1d_invprob_lin(DataFile,UseRecursive,ExaCoefFile)
% input: 
%     DataFile: text file with the potential drop data and frequencies
%     UseRecursive: index (=1 if we use recursive procedure, = 0 if not)
%     ExaCoefFile: optional, exact coefficient profiles for simulated data only        


ParameterFile = 'linear_invprob_parameters.txt';
HomParFile = 'hom_invprob_parameters.txt';

OutputFile = [DataFile(1:end-4),'_result.txt'];

DataFile % show the data file name

[rho11,rho12,rho21,rho22,PerturbationParameter,Nz,Depth,RegPar] = acpd1d_invprob_lin_get_parameters(ParameterFile);



if UseRecursive   % run the linearized inverse problem with recursive procedure: 
    figname1 = [DataFile(1:end-4),'_rec_cond',num2str(Nz),'.eps'];
    figname2 = [DataFile(1:end-4),'_rec_prod',num2str(Nz),'.eps'];

    [Conductivity,ProductSigmaMu] = acpd1d_invprob_lin_recursive(ParameterFile,HomParFile,DataFile);

else    % run the linear inverse problem without recursive procedure:
    figname1 = [DataFile(1:end-4),'_cond',num2str(Nz),'.eps'];
    figname2 = [DataFile(1:end-4),'_prod',num2str(Nz),'.eps'];

    [Freq,PotentialDrop] = acpd1d_invprob_load_data(DataFile); 
    [Sigma0,Mu0] = acpd1d_invprob_hom(HomParFile,Freq,PotentialDrop);
    Sigma0
    Mu0
    [Conductivity,ProductSigmaMu] = acpd1d_invprob_lin(ParameterFile,DataFile,Sigma0,Mu0);
    % ProductSigmaMu: the product of the conductivity sigma and permeability mu.
end


Permeability = ProductSigmaMu/Conductivity; 

%save the result to a MAT file for later use: 
Result = [Conductivity.Depth, Conductivity.CoefValue, ProductSigmaMu.CoefValue];
dlmwrite(OutputFile,Result,'delimiter',' ');

% compare the simulated data with the measured data:
acpd1d_invprob_compare_result_data(Conductivity,Permeability,Freq,PotentialDrop,rho11,rho12,rho21,rho22,DataFile);

    
if nargin < 2
    acpd1d_plot_coef(Depth,Conductivity,ProductSigmaMu,figname1,figname2);
else 
    acpd1d_plot_coef(Depth,Conductivity,ProductSigmaMu,figname1,figname2,ExaCoefFile);
end