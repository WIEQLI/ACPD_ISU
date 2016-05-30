function [Sigma,Mu] = acpd1d_invprob_hom(ParameterFile,Freq,PotentialDrop,SigmaInitial,MuInitial)
% solve the linearized inverse problem for the AC potential drop. 
% Input: ParameterFile: text file containing input parameters
%        Freq: list of frequencies used
%        PotentialDrop: potential drop data
%     Output: Conductivity and Relative Permeability, two values
% =========================================================================

% options for the MATLAB optimization algorithm
options = optimset('MaxIter',100,'TolFun',1e-15);
options = optimset(options,'GradObj','off','display','iter','Algorithm','sqp','TolX',1e-14,'MaxFunEvals',1000,'MaxIter',200);


MU0 = 12.5663706143592e-07; % permeabliity of the free space

% ---- load input parameters: 
fid = fopen(ParameterFile);
if fid == -1
    error('The parameter file not found');
end
rho11 = fscanf(fid,'%f',1); rho12 = fscanf(fid,'%f',1); 
rho21 = fscanf(fid,'%f',1); rho22 = fscanf(fid,'%f',1); fgetl(fid);
SigmaInit = fscanf(fid,'%f',1); SigmaLower = fscanf(fid,'%f',1); SigmaUpper = fscanf(fid,'%f',1); fgetl(fid);
MuInit = fscanf(fid,'%f',1); MuLower = fscanf(fid,'%f',1); MuUpper = fscanf(fid,'%f',1); fgetl(fid);
fclose(fid); % end of loading the parameters:

% if the initial guess for Sigma and Mu are provided in the code, use them:
if nargin > 3
    SigmaInit = SigmaInitial;
    MuInit = MuInitial; 
end

% solving the optimization problem:

Mu0 = (MuLower + MuUpper)/2*MU0; % average of lower and upper bounds as the scaling factor
Sigma0  = (SigmaLower + SigmaUpper)/2;
Scaling = [Sigma0; Mu0];

InitGuess = [SigmaInit; MuInit*MU0]./Scaling; 
lb = [SigmaLower; MuLower*MU0]; % lower bound
ub = [SigmaUpper; MuUpper*MU0]; % upper bound

lb = lb./Scaling;
ub = ub./Scaling;
Sol = fmincon(@(V)objfun(V),InitGuess,[],[],[],[],lb,ub,[],options);
% Sol = fminunc(@(V)objfun(V),InitGuess,options);

Sigma = Sol(1,end)*Sigma0; Mu = Sol(2,end)*Mu0;
Mu = Mu/MU0; % relative permeability

% objective function: 
    function F = objfun(V)
        sigma0 = V(1)*Sigma0; mu0 = V(2)*Mu0; 
        PD = acpd1d_function_D0(Freq,sigma0,mu0,rho11, rho12, rho21,rho22); 
        F = sum((real(PD) - real(PotentialDrop)).^2)  +  sum((imag(PD) - imag(PotentialDrop)).^2);
        F = F*10^6;         
    end
end