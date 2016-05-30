function D0 = acpd1d_function_D0(freq,sigma0,mu0,rho11, rho12, rho21,rho22)
% calculate the zero order term in the potential drop. This term
% corresponds to the potential drop of a homogeneous medium
% freq: a vector of frequencies

% The formula used here gives the opposite sign as the one in the written
% file. 

if (abs(rho11 - rho22) < 2*eps) && (abs(rho12 - rho21) < 2*eps)
    Symmetric = true;
else
    Symmetric = false;
end

% calculate D0 through the explicit formula of the potential drop:
Nf = length(freq);

D0 = zeros(Nf,1);
 
Tolerance = 1e-6*abs(1/rho11 - 1/rho12);

for n = 1:Nf

    k = sqrt(1i*2*pi*freq(n)*sigma0*mu0); 
    integrand = @(x)exp(1i*k*x)./x.^2;
    
    if Symmetric
        D0(n) = 1/pi/sigma0*(integral(integrand,rho11,rho12,Tolerance) - 1i*k*log(rho12/rho11));        
    else
        D0(n) = 1/2/pi/sigma0*(integral(integrand,rho11,rho12,Tolerance) - 1i*k*log(rho12/rho11))...
              + 1/2/pi/sigma0*(integral(integrand,rho22,rho21,Tolerance) - 1i*k*log(rho21/rho22));
    end
end


% % Calculate D0 through the integral from 0 to infinity:
% k = sqrt(2*pi*freq*sigma0*mu0); 
% 
% dkappa = 0.1; alpha = 10; % these parameters are chosen by trial-and-error.
% kappa = 0.000001:dkappa:alpha*k;
% 
% Kernel = besselj(0,kappa.*rho11)-besselj(0,kappa.*rho12)-besselj(0,kappa.*rho21)+besselj(0,kappa.*rho22);
% 
% I = dkappa*sum( (sqrt(1 - 1i*k^2./kappa.^2) - 1).*Kernel) + 1/rho11 + 1/rho22 - 1/rho21 - 1/rho12;
% D02 = I/2/pi/sigma0
% 
% 



