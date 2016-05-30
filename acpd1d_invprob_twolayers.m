function Sol = acpd1d_invprob_twolayers(meadat, Freq, rho11,rho12,rho21,rho22,InitGuess)


Nf = length(Freq);

options = optimset('MaxIter',20,'TolFun',1e-9);
options = optimset(options,'GradObj','off','display','iter','Algorithm','sqp');
RegPar = 1e-11;

Sigma0 = InitGuess(1); 
Sigma1 = InitGuess(2);
Mu0 = InitGuess(3);
Mu1 = InitGuess(4);
Beta = InitGuess(5);

InitGuess = ones(5,1);
lb = zeros(5,1); ub = 10*ones(5,1);

if nargin > 7 % constrained optimization  
    Sol = fmincon(@(V)objfun(V),InitGuess,[],[],[],[],lb,ub,[],options);
else % unconstrained problem
    Sol = fminunc(@(V)objfun(V),InitGuess,options);
end

Sol = Sol(:,end);
Sol(1) = Sol(1)*Sigma0;
Sol(2) = Sol(2)*Sigma1;
Sol(3) = Sol(3)*Mu0;
Sol(4) = Sol(4)*Mu1;
Sol(5) = Sol(5)*Beta;


    function F = objfun(V)
        Sig0 = Sigma0*V(1);
        Sig1 = Sigma1*V(2);
        MU0 = Mu0*V(3);
        MU1 = Mu1*V(4);
        Zmin = Beta*V(5);
        
        F = 0; 
        for n = 1:Nf
            PD = acpd1d_pd_twolayers(Zmin, Freq(n), MU0,MU1,Sig0,Sig1,rho11,rho12,rho21,rho22);
            F = F + (real(PD) - real(meadat(n))).^2 ;% +  (imag(PD) - imag(meadat(n))).^2;
        end
        F = F + RegPar*sum((V - InitGuess).^2);
       F = F*10000;
    end
end