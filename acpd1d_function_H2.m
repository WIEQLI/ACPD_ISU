function H2 = acpd1d_function_H2(freq,sigma0,mu0,rho11, rho12, rho21,rho22,y)
% calculate the value of the integral H2 in the AC potential drop linearized
% formulation. 
% H2 = 1/2/pi/sigma0*int_{0}^infty exp() Omega(kappa)/kappa dkappa
% where Omega(kappa) = J_0(rho11*kappa) - J_0(rho12*kappa) - J_0(rho21*kappa) + J_0(rho22*kappa)
% J_0: Bessel function of order zero.
% Each column of H2 is the value at a fixed frequency
    
% Numerical approximation    

if size(y,2) > 1
    y = y.'; % y must be a column vector. 
    row = 1;
end

omega = 1i*2*pi*freq*sigma0*mu0; 
dkappa = 10; alpha = 100000; 
kappa = 0.0000001:dkappa:alpha;

Kernel = besselj(0,kappa.*rho11)-besselj(0,kappa.*rho12)-besselj(0,kappa.*rho21)+besselj(0,kappa.*rho22);

Nf = length(freq);
Ny = length(y);
Nk = length(kappa);

y = y*ones(1,Nk);
kappa = ones(Ny,1)*kappa;
Kernel = ones(Ny,1)*Kernel;

I22 = zeros(Ny,Nf);
for n = 1:Nf  
    ksq = omega(n);
    integrand = exp(sqrt(kappa.^2 - ksq).*y*2)./kappa.*Kernel;
%     figure(1); plot(kappa,real(integrand)*freq^2); title('real part'); 
%     figure(2); plot(kappa,imag(integrand)*freq^2); title('imaginary part');
    I22(:,n) = sum(integrand,2)*dkappa;
end

H2 = I22*ksq/2/pi/sigma0;
if row==1
    H2  = H2.';
end

