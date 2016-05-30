% test the numerical simulation of the potential drop computation: 
% verify the validity of computation of the potential drop for a two layered medium:

MU0 = 4*pi*10^(-7); % permeability of vacuum 

% % data set 1: 
% dataset = 'ACPD1016surf'; 
% data = dlmread([dataset,'.txt']); % data file from Rashid
% sigma0 = 4.42*10^6;  % second layer extended to -infinity
% sigma1 = 2.865*10^6; % surface layer
% mu0 = 67.046*MU0; % permeability of the second layer
% mu1 = 31.821*MU0; % permeability of the surface layer
% Zmin = -0.0014; % the depth of the surface layer in meters

% % data set 2: 
% dataset = 'ACPD8002surf'; 
% data = dlmread([dataset,'.txt']); % data file from Rashid
% sigma0 = 3.428*10^6;  % second layer extended to -infinity
% sigma1 = 2.115*10^6; % surface layer
% mu0 = 62*MU0; % permeability of the second layer
% mu1 = 22.443*MU0; % permeability of the surface layer
% Zmin = -0.0009; % the depth of the surface layer in meters

% data set 3: 
dataset = 'ACPD1p5mm176M1'; 
data = dlmread([dataset,'.txt']); % data file from Rashid
sigma0 = 2.85*10^6; 
sigma1 = 2.663*10^6; % surface layer 2.663
mu0 = 75.64*MU0; % permeability of the second layer
mu1 = 48.087*MU0; % permeability of the surface layer 48.087
Zmin = -0.0018; % the depth of the surface layer in meters
Permittivity = PwConstCoefficient();

% frequencies and measured potential drop:
FreIdx = size(data,1)-5:size(data,1); % choose a set of frequencies
PD_mea = data(FreIdx,8) + 1i*data(FreIdx,9);
Freq = data(FreIdx,1); 


rho11 = 0.0015; rho12 = 0.003; rho21 = 0.003; rho22 = 0.0015; 

% % Estimate the parameters: 
% dataset = [dataset,'_2'];
% InitGuess = [sigma0; sigma1; mu0; mu1; Zmin];
% Sol = acpd1d_invprob_twolayers(PD_mea, Freq, rho11,rho12,rho21,rho22,InitGuess);
% 
% sigma0 = Sol(1)
% sigma1 = Sol(2)
% mu0 = Sol(3); disp(mu0/MU0)
% mu1 = Sol(4); disp(mu1/MU0)
% Zmin = Sol(5)

% Estimate the surface layer: 

InitGuess = [sigma1; mu1];
lb = [1e6; MU0]; ub = [1e7; 150*MU0];

[Sigma,Mu] = acpd1d_invprob_hom(PD_mea, Freq, rho11,rho12,rho21,rho22,InitGuess)

RelativeMu = Mu/4/pi*10^7

sigma0 = Sigma
mu0 = Mu;



% Simulate the data: 
z = linspace(Zmin,0,2); % to aproximate the integral with kernel I2 by a finite sum
EpsSigma = (sigma1/sigma0 - 1) + z(1:end-1)*0; %epsilon*alpha(z)
Alpha = EpsSigma/(1 + EpsSigma(1)); % the unknown coefficient in the linearized problem
EpsMu = (mu1/mu0 - 1) + z(1:end-1)*0;
Eta = EpsMu/(1 + EpsSigma(1)); %unknown of permeability

% Freq = (10000)';
PD_lin = 0*Freq; 
% PD_mea = PD_lin;





figure; set(gca,'fontsize',15); 
plot(Freq,real(PD_lin),'-b','linewidth',2); 
hold on; plot(Freq,real(PD_mea),'--r','linewidth',2); hold off;
xlabel('Frequency (Hz)'); ylabel('Re \Delta v (V)'); 
legend('Simulation','Measured data');
title('Real part of the potential drop');
print('-depsc2',['2layers_Re_',dataset,'.eps']);


figure; set(gca,'fontsize',15); 
plot(Freq,imag(PD_lin),'-b','linewidth',2); 
hold on; plot(Freq,imag(PD_mea),'--r','linewidth',2); hold off;
xlabel('Frequency (Hz)');  ylabel('Im \Delta v (V)'); 
legend('Simulation','Measured data');
title('Imaginary part of the potential drop');
print('-depsc2',['2layers_Im_',dataset,'.eps']);

