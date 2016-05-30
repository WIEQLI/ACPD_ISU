function acpd1d_invprob_compare_result_data(Conductivity,Permeability,Freq,PotentialDrop,rho11,rho12,rho21,rho22,DataFileName)
% compare the simulated data with the estimated coefficients, Conductivity
% and Permeability, and the measured data. 
% this function is called in the inverse routines. 



PD_sim =  acpd1d_potentialdrop(Conductivity,Permeability,Freq,rho11,rho12,rho21,rho22,'linear'); % measured data

Freq = log(Freq);

figure(3); set(gca,'fontsize',18); 
plot(Freq,real(PD_sim),'-b','linewidth',2);
hold on; plot(Freq,real(PotentialDrop),'--r','linewidth',2);
hold off;
title('Real part of the potential drop'); 
legend('Simulation','Measured data');
xlabel('Frequency in log scale'); ylabel('Re(\Delta v)');
print('-depsc2',[DataFileName(1:end-4),'_pd_re.eps']);

figure(4); set(gca,'fontsize',18); 
plot(Freq,imag(PD_sim),'-b','linewidth',2);
hold on; plot(Freq,imag(PotentialDrop),'--r','linewidth',2);
hold off;
title('Imaginary part of the potential drop'); 
legend('Simulation','Measured data');
xlabel('Frequency in log scale'); ylabel('Img(\Delta v)');
print('-depsc2',[DataFileName(1:end-4),'_pd_im.eps']);

