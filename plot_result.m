% Plot the result:



% [ExSigma,ExMu] = acpd1d_simulate_data('simulation_parameters.txt');
[Sigma0,Mu0] = acpd1d_invprob_hom('hom_invprob_parameters.txt',DataFile);
[EstSigma,EstProd] = acpd1d_invprob_lin('linear_invprob_parameters.txt',DataFile,OutputFile,Sigma0,Mu0);
figname1 = 'conductivity_15.eps';
figname2 = 'prod_cond_perm_15.eps';


ExProd = ExSigma*ExMu; % product of sigma and mu.
Zmin = -0.003; 

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

