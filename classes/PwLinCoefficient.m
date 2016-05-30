% a subclass of the Coefficient class: 
% PwLinCoefficient is a class which handles piecewise linear coefficients.
% Differences from the superclass: the functions "plot" and
% "WeightedIntegral"
% updated March 25, 2016. 

classdef PwLinCoefficient < Coefficient    
    
    methods
        function coef = PwLinCoefficient(depth,coefval,refval)
            coef@Coefficient(depth,coefval,refval);
            
        end      
        
        function plot(Coef,fig,colorcode)
            if nargin > 1
                figure(fig); 
            else
                figure;
                set(gca,'fontsize',15);
            end
            if nargin < 3
                colorcode = '-b'; % blue
            end
            set(gca,'fontsize',18);
            plot([Coef.Depth(1)-abs(Coef.Depth(2)-Coef.Depth(1)), Coef.Depth(1)],[Coef.CoefValue(1), Coef.CoefValue(1)],colorcode,'linewidth',2);
            hold on;
            plot(Coef.Depth,Coef.CoefValue,colorcode,'linewidth',2);
            xlabel('z'); 
            hold off;
        end
        
        function newcoef = Interpolation(Coef,z) % linear interpolation. z can be a vector
            newcoef = Coef;
            newcoef.Depth = z;
            if z(1) < Coef.Depth(1) % extrapolation
                Coef.Depth = [z(1); Coef.Depth];
                Coef.CoefValue = [Coef.CoefValue(1); Coef.CoefValue];
            end
            newcoef.CoefValue = interp1(Coef.Depth,Coef.CoefValue,z,'linear');
            newcoef.RelativeCoef = newcoef.CoefValue/newcoef.RefCoef;            
        end
        
        function CoefVec = CoefficientVector(Coef,Weight,Zmin) % calculate the integral of f*Coef from Depth(idx) to Depth(idx+1)
            
            tolerance = 1e-8; % tolerance for numerical integration.
            
            Nz = length(Coef.Depth); % number of unknowns
            CoefVec  = zeros(1,Nz); % a: mean of H1
            Z = Coef.Depth;
            
            CoefVec(1) = integral(Weight,Zmin,Z(1),tolerance) ...
                       + integral(@(y)(Weight(y).*(Z(2)-y)),Z(1),Z(2),tolerance)/(Z(2) - Z(1));
            for n = 2:Nz-1
                 CoefVec(n) = integral(@(y)(Weight(y).*(y-Z(n-1))),Z(n-1),Z(n),tolerance)/(Z(n)-Z(n-1)) ...
                            + integral(@(y)(Weight(y).*(Z(n+1)-y)),Z(n),Z(n+1),tolerance)/(Z(n+1)-Z(n));
            end
            CoefVec(Nz) = integral(@(y)(Weight(y).*(y-Z(Nz-1))),Z(Nz-1),Z(Nz),tolerance)/(Z(Nz)-Z(Nz-1));
        end
        
    end
        
end
