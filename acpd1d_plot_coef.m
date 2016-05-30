function acpd1d_plot_coef(Depth,EstSigma,EstProd,figname1,figname2,ExaCoefFile)
% plot the conductivity and product of conductivity and permeability profiles in
% depth for the inverse problems. 
% compare with the exact coefficient values for simulated data.

Zmin = Depth(1)*1.2; 

if nargin > 5 && exist(ExaCoefFile,'file') % for simulated data with avaiable exact coefficients:
    eval(['load ' ExaCoefFile]); % load the exact coefficients

    ExProd = ExSigma*ExMu; % product of sigma and mu.

    ExSigma.plot_coef(1,'--r',Zmin);
    hold on; 
    EstSigma.plot_coef(1,'-b',Zmin); hold off;

    legend('Exact coefficient','Estimated coefficient');
    xlabel('z (m)'); ylabel('\sigma(z)');
    axis([Zmin 0 min(min(ExSigma.CoefValue),min(EstSigma.CoefValue))*0.9, max(max(ExSigma.CoefValue),max(EstSigma.CoefValue))*1.1]);
    print('-depsc2',figname1);

    ExProd.plot_coef(2,'--r',Zmin); hold on;
    EstProd.plot_coef(2,'-b',Zmin); hold off;
    legend('Exact coefficient','Estimated coefficient');
    xlabel('z (m)'); ylabel('\sigma(z)\mu(z)');
    axis([Zmin 0 min(min(ExProd.CoefValue),min(EstProd.CoefValue))*0.9, max(max(ExProd.CoefValue),max(EstProd.CoefValue))*1.1]);
    print('-depsc2',figname2);

else 

    EstSigma.plot_coef(1,'-b',Zmin); hold off;

    xlabel('z (m)'); ylabel('\sigma(z)');
    axis([Zmin 0 min(EstSigma.CoefValue)*0.9, max(EstSigma.CoefValue)*1.1]);
    print('-depsc2',figname1);

    EstProd.plot_coef(2,'-b',Zmin); hold off;
    xlabel('z (m)'); ylabel('\sigma(z)\mu(z)');
    axis([Zmin 0 min(EstProd.CoefValue)*0.9, max(EstProd.CoefValue)*1.1]);
    print('-depsc2',figname2);
end