function pd = acpd1d_potentialdrop(Conductivity,Permeability, Freq, rho11,rho12,rho21,rho22,method)
% calculate the potential drop
% Conductivity, Permeability: two objects of class "Coefficient"
% Freq: set of frequencies
% rho11,rho12,rho21,rho22: distances between pins in the experiment
% method: could be "IVP" for using the IVP solver, "twolayer" for two
% layered model, "hom" for homogeneous medium, "linear" for the linearized
% model. Default method: "IVP".


if nargin < 8
    method = 'IVP';
end
if max(abs(Permeability.Depth - Conductivity.Depth)) > eps
    error('Conductivity and permeability are not consistent');
end


if strcmpi(method,'IVP')
    disp('Default method used: IVP');  
    pd = potentialdrop_IVP(Conductivity,Permeability,Freq,rho11,rho12,rho21,rho22);
    
elseif strcmpi(method,'hom') && (length(Conductivity.Depth) == 1) && (length(Permeability.Depth) == 1)
    pd = acpd1d_function_D0(Freq,Conductivity.CoefValue,Permeability.CoefValue,rho11, rho12, rho21,rho22); 
    
elseif strcmpi(method,'twolayer') && (length(Conductivity.Depth) == 2) && (length(Permeability.Depth) == 2)
    pd = potentialdrop_twolayer(Conductivity,Permeability,Freq,rho11,rho12,rho21,rho22);

elseif strcmpi(method,'linear') 
    pd = potentialdrop_linear(Conductivity,Permeability,Freq,rho11,rho12,rho21,rho22);
else
    disp('The requested method is not supported. Choose between hom, twolayer, linear, and IVP.');    
end


    % subfunctions: two-layer model: 
    function PD = potentialdrop_twolayer(Conductivity,Permeability,Freq,rho11,rho12,rho21,rho22)    
        NumFreq = length(Freq); 
        Zmin = Conductivity.Depth(1);
        PD = zeros(NumFreq,1);
        [sigma0,sigma1] = Permivivity.CoefValue;
        [mu0,mu1] = Permeability.CoefValue;

        dkappa = 0.5; Alpha = 5000; % these parameters are chosen by trial-and-error.
        kappa = 0.0000001:dkappa:Alpha;
        Kernel = besselj(0,kappa.*rho11)-besselj(0,kappa.*rho12)-besselj(0,kappa.*rho21)+besselj(0,kappa.*rho22);

        for n = 1:NumFreq  
            freq = Freq(n);
            
            ksq = 1i*2*pi*freq*sigma0*mu0; 

            gamma = sqrt(1- ksq./kappa.^2);
            gamma1 = sqrt((1- ksq*sigma1/sigma0*mu1/mu0./kappa.^2));
            v1 = gamma1./gamma; 
            v2 = sigma1/sigma0;
            exponent = exp(2*gamma1.*kappa*Zmin);
            D = (v1 + v2 + (v2 - v1).*exponent)./(v1 + v2 + (v1 - v2).*exponent);
            integrand = (D-1).*gamma1.*Kernel;

            PD(n) = dkappa*sum(integrand)/2/pi/sigma1 + acpd1d_function_D0(freq,sigma1,mu1,rho11, rho12, rho21,rho22);
        end
    end

    % subfunction: linear model: 
    function PD = potentialdrop_linear(Conductivity,Permeability,Freq,rho11,rho12,rho21,rho22)  
        
        sigma0 = Conductivity.RefCoef;
        mu0 = Permeability.RefCoef;
        
        D0 = acpd1d_function_D0(Freq,sigma0,mu0,rho11, rho12, rho21,rho22);
        D0 = D0(:);

        [A,B] = CoefficientMatrix(Conductivity,Freq,sigma0,mu0,rho11,rho12,rho21,rho22); 
        % updated: A and B correspond to two coefficients: mu*sigma and sigma)
        
        Alpha = (Conductivity.RelativeCoef-1)/(Conductivity.RelativeCoef(end));
        P = (Conductivity.CoefValue.*Permeability.CoefValue); %product of mu and sigma (we combine the two in one coefficient)
%        Eta = (Permeability.RelativeCoef-1)/(Conductivity.RelativeCoef(end));
        Eta = (P/(sigma0*mu0) -1)/(Conductivity.RelativeCoef(end)); %Eta = (mu*sigma/(mu0*sigma0) - 1)/(1 + epsilon*alpha(0));
                
        PD = D0 + A*Alpha + B*Eta;         
    end

    function PD = potentialdrop_IVP(Conductivity,Permeability,Freq,rho11,rho12,rho21,rho22)   
        
        % NEED UPDATE
        NumFreq = length(Freq); 
        Zmin = Conductivity.Depth(1);
        PD = zeros(NumFreq,1);
        sigma0 = Conductivity.CoefValue(1); % value at the deep depth
        mu0 = Permeability.CoefValue(1);
        sigma0mu0 = sigma0*mu0;

        dkappa = 5; Alpha = 40000; % these parameters are chosen by trial-and-error.
        kappa = 0.0000001:dkappa:Alpha; % range of kappa approximating (0, infinity)
        Kernel = besselj(0,kappa.*rho11)-besselj(0,kappa.*rho12)-besselj(0,kappa.*rho21)+besselj(0,kappa.*rho22);

        % discretization parameters: mesh
        dz = 1/sqrt(Alpha^2 + max(Freq)*2*pi*sigma0mu0)/2; % step size in depth, choose small for stability of the numerical scheme.
        dzsq = dz^2; % square of dz
        z = Zmin:dz:0; Nz = length(z); 
        % correct the mesh so that it contains 0 as the end point. 
        Nz = Nz + 1; z = linspace(Zmin,0,Nz); dz = z(2)-z(1);
        
        % Linaerly interpolate the coefficients: 
        CondInt = Conductivity.Interpolation(z);
        PermInt = Permeability.Interpolation(z); 
        
        MuSigma = CondInt.CoefValue.*PermInt.CoefValue; % product of sigma and mu.
        p = MuSigma/sigma0mu0; % p(z) = mu(z)*sigma(z)/mu0/sigma0
        logs = log(CondInt.CoefValue/sigma0); % ln(s), where s = sigma/sigma0.
        
       
        NumKappa = length(kappa); 
        
        for n = 1:NumFreq
            freq = Freq(n);
            DerF_Zero = zeros(1,NumKappa);       
            for nk = 1:NumKappa
                kapp = kappa(nk);
                
                ksq = 1i*2*pi*freq*sigma0mu0; 
                gamma0 = sqrt(kapp^2 - ksq);  
                gammasq = kapp^2 - ksq*p; 
                
                F = zeros(1,Nz); % the function F
                Fder = zeros(1,Nz); % derivative of F. 
                
                % initial condition at z = Zmin:
                F(1) = exp(Zmin*gamma0); Fder(1) = gamma0*exp(Zmin*gamma0);                
                % iteration: 
                for nz = 2:Nz                 
                   A = 1 - dzsq*gammasq(nz)/4 - (logs(nz)-logs(nz-1))/2; 
                   B = 2 - A;
                   C = dz/2*(gammasq(nz) + gammasq(nz-1));                    
                   Fder(nz) = (B*Fder(nz-1) + C*F(nz-1))/A;
                   F(nz) = F(nz-1) + (Fder(nz) + Fder(nz-1))*dz/2;
                end                
                C = 1/2/pi/F(Nz); % note that this C is not the same as in the for loop. 
                DerF_Zero(nk) = Fder(Nz)*C; %F'(0,omega,kappa). 
            end            
            PD(n) = -sum(Kernel.*DerF_Zero./kappa)*dkappa/CondInt.CoefValue(Nz); %9/30/2015: accuracy must be checked. The integrand's decaying is not fast enough            
        end
    end    
end






