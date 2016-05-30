function H1 = acpd1d_function_H1(freq,sigma0,mu0,rho11, rho12, rho21,rho22,y)
% calculate the value of the integral H1 in the AC potential drop linearized
% formulation. 
% H1 = 1/pi/sigma0*int_{0}^infty kappa exp(..)Omega(kappa) dkappa
% where Omega(kappa) = J_0(rho11*kappa) - J_0(rho12*kappa) - J_0(rho21*kappa) + J_0(rho22*kappa)
% J_0: Bessel function of order zero.
% Output is a column vector or a matrix, each column is the function H1 at each frequency.    
% Explicit formula.    
    

if size(y,2) > 1
    y = y.'; % y must be a column vector. 
    row = 1;
end

k = sqrt(2*pi*freq*sigma0*mu0); 

% the first part (the non-decaying part):
H1 = hankel_transform_H(-2*y,k,rho11) - hankel_transform_H(-2*y,k,rho21) ...
    - hankel_transform_H(-2*y,k,rho12) + hankel_transform_H(-2*y,k,rho22);

H1 = H1/pi/sigma0;
if row==1
    H1 = H1.';
end

    % subfunction:
    function F = hankel_transform_H(a, k, rho)
    % calculate the expression:
    % a*(rho^2 + a^2)^(-3/2) * exp(-k*exp(-1i*pi/4)*(rho^2 + a^2)^(1/2)) *(1 + k*exp(-1i*pi/4)*(rho^2 + a^2)^(1/2))
    % This expression should be equal to the Hankel transform: 
    %      int_0^infty kappa*exp(-a*(kappa^2 - 1ik^2)^(1/2)) + J_0(kappa*rho) dkappa. 
    % a: a column vector, k: a row vector
    % F: a matrix, each column is at a frequency.
    %%%%

    if size(a,2) > 1
        a = a';
    end
    if size(k,1) > 1
        k = k';
    end

    Na = length(a); 
    Nk = length(k);

    a = a*ones(1,Nk); % convert to matrix
    k = ones(Na,1)*k;
    F = a.*((rho^2 + a.^2).^(-3/2)).*exp(-k.*exp(-1i*pi/4).*(rho^2 + a.^2).^(1/2)).*(1 + k.*exp(-1i*pi/4).*(rho^2 + a.^2).^(1/2));
    end
end

