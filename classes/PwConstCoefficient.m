% a subclass of the coefficient class, used to handle piecewise constant
% coefficients

classdef PwConstCoefficient < Coefficient
    
    methods
        function coef = PwConstCoefficient(depth,coefval,refval)
            coef@Coefficient(depth,coefval,refval);
            
        end
      
        function plot_coef(Coef,fig,colorcode,Zmin)
            
            N = length(Coef.Depth);            
            if nargin > 1
                figure(fig); 
            else
                figure;
                set(gca,'fontsize',15);
            end
            if nargin < 3
                colorcode = '-b'; % blue
            end
            if nargin < 4
                Zmin = Coef.Depth(1)-abs(Coef.Depth(2)-Coef.Depth(1));
            end
            plot([Zmin, Coef.Depth(1)],[Coef.CoefValue(1), Coef.CoefValue(1)],colorcode,'linewidth',2);
            hold on;
            for n = 2:N
                plot([Coef.Depth(n-1), Coef.Depth(n)],[Coef.CoefValue(n),Coef.CoefValue(n)],colorcode,'linewidth',2);
            end
            hold off;
            
            xlabel('z'); 
        end
        
        function CoefVec = CoefficientVector(Coef,Weight,Zmin)
        
            tolerance = 1e-8;
            
            Nz = length(Coef.Depth); % number of unknowns
            CoefVec  = zeros(1,Nz); % a: mean of H1
            Z = Coef.Depth;
            CoefVec(1) = integral(Weight,Zmin,Z(1),tolerance);
            for n = 2:Nz
                 CoefVec(n) = integral(Weight,Z(n),Z(n+1),tolerance);
            end
        
        end
    end   
    
end
