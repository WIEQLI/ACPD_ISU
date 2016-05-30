% define class "Coefficient" for coefficient of the 1d AC potential drop
% problem. 
% Properties: 
%   - Depth: a column vector of depth with negative values, starting from smallest value to
%   zero. 
%   - CoefValue: Value of the coefficient associated with the depth. The
%   relationship between them depends on the type of coefficient being used,
%   defined in subclasses. The number of elements of CoefValue and Depth are
%   equal. Coefficient(z) = Coefficient(1) for z < Depth(1)
%   - RefCoef: reference value (sigma_0 or mu_0). 
%   - RelativeCoef = CoefValue/RefCoef
% updated Oct 6 2015: add multiplication method mtimes: overloading of the
% element-by-element multiplication operator, *. 
% Updated March 25, 2016: add more comments for easier to read. 




classdef Coefficient
    properties (Access = public)
        Depth         % a vector of depths. 
        CoefValue     % absolute coefficient value
        RefCoef % a constant reference value at which the perturbation of the coefficient is taken.
        RelativeCoef % relative coefficient value = CoefValue/RefCoef. 
    end
 
    methods (Access = public)
        function Coef = Coefficient(depth,val,refcoef) % construction
            if nargin == 3
                if isnumeric(depth) && isnumeric(val) && (length(depth) == length(val))
                    depth = depth(:); % convert to a column vector
                    val = val(:); % convert to a column vector
                    Coef.Depth = depth;
                    Coef.CoefValue = val;
                    Coef.RelativeCoef = val/refcoef;
                    Coef.RefCoef = refcoef;
                    
                else
                    disp('Error in Coefficient class construction: The input arguments are not valid');
                end
                
            end
        end
        
                
        function plot_coef(Coef,fig,colorcode,Zmin) 
            if nargin > 1
                figure(fig); 
            else
                figure;
                set(gca,'fontsize',15);
            end
            if nargin < 3
                colorcode = '*b'; % blue
            end
            if nargin < 4
                Zmin = Coef.Depth(1)-abs(Coef.Depth(2)-Coef.Depth(1)); % if Zmin is not provided, just add one more mesh size. 
            end
            set(gca,'fontsize',18);

            Coef.Depth = [Zmin; Coef.Depth];
            Coef.CoefValue = [Coef.CoefValue(1); Coef.CoefValue];
            plot(Coef.Depth,Coef.CoefValue,colorcode,'linewidth',2);
            xlabel('z'); 
        end
        
        newCoef = Interpolation(Coef,z) % interpolate the value of Coefficient to a grid given in vector "z". 
        % The implementation of this function depends on the type of coefficient (piecewise linear, piecewise constant, ...). 
        % Therefore its implementation is done in subclasses. 

        function Coef3 = mtimes(Coef,Coef2) % element-by-element multiplication (overload operator *)
            if length(Coef.Depth) ~= length(Coef2.Depth)
                error('Coefficient multiplication error: their sizes are not equal');
            end
            Coef3 = Coef; 
            Coef3.CoefValue = Coef.CoefValue.*Coef2.CoefValue;
            Coef3.RefCoef = Coef.RefCoef*Coef2.RefCoef; 
            Coef3.RelativeCoef = Coef.RelativeCoef.*Coef2.RelativeCoef; 
        end
            
        function Coef3 = mrdivide(Coef,Coef2)  % element-by-element division (overload the operator /)
            if length(Coef.Depth) ~= length(Coef2.Depth)
                error('Coefficient multiplication error: their sizes are not equal');
            end
            if min(abs(Coef.CoefValue)) < eps
                error('Division by zero!');
            end
            Coef3 = Coef; 
            Coef3.CoefValue = Coef.CoefValue./Coef2.CoefValue;
            Coef3.RefCoef = Coef.RefCoef/Coef2.RefCoef; 
            Coef3.RelativeCoef = Coef.RelativeCoef./Coef2.RelativeCoef; 
        end
  
        CoefVec = CoefficientVector(Coef,Weight,Zmin)  % calculate the integral of f*Coef in each subinterval
        % how this function is implemented depends on the type of the coefficient. Therefore its implementation is done
        % in subclasses. 
        

        function [A,B] = CoefficientMatrix(Coef,Freq,sigma0,mu0,rho11,rho12,rho21,rho22)
            % calculate the coefficient matrices of the linearized model.
            % The linearized model reads: A\alpha + B\eta = RHS. 
            % eta corresponds to the product of Mu and Sigma
            
            
            Nz = length(Coef.Depth); % number of unknowns
            Nf = length(Freq); % number of frequencies
            A = zeros(Nf,Nz); % coefficient matrix w.r.t. conductivity
            B = zeros(Nf,Nz); % coefficient matrix w.r.t. permeability

            for m = 1:Nf
                freq = Freq(m);
                Zmin = -10*sqrt(1/pi/freq/mu0/sigma0); % 10 times of skin depth. The kernel is assumed to be negligible at this depth
                
                H1 = @(y)acpd1d_function_H1(freq,sigma0,mu0,rho11, rho12, rho21,rho22, y);
                H2 = @(y)acpd1d_function_H2(freq,sigma0,mu0,rho11, rho12, rho21,rho22, y);
                CoefVec1 = CoefficientVector(Coef,H1,Zmin);
                CoefVec2 = CoefficientVector(Coef,H2,Zmin);
                
                A(m,:) = CoefVec1 - 2*CoefVec2; % this matrix corresponds to the linearization of p(z) = mu*sigma.
                B(m,:) = CoefVec2;
            end            
        end
        
        
    end % end of methods    
    
end  % end of class
    
    
    
    
    
    