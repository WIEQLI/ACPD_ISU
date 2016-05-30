function [ExSigma,ExMu] = acpd1d_simulate_data(parameterfile,OutputCoefFile)
% simulate the data for the acpd problem. 
% parameterfile: the file containing the parameters of the simulation,
% OutputCoefFile: a .mat file containing the coefficients classes (for
% visualization and comparison with the reconstructed coefficients

MU0 = 12.5663706143592e-07; % permeabliity of the free space

% ---- load simulation parameters: 
fid = fopen(parameterfile);
if fid == -1
    error('The parameter file not found');
end
OutputFileName = fscanf(fid,'%s',1); fgetl(fid);   

rho11 = fscanf(fid,'%f',1); rho12 = fscanf(fid,'%f',1); 
rho21 = fscanf(fid,'%f',1); rho22 = fscanf(fid,'%f',1); 
fgetl(fid);

NumGridPoints = round(fscanf(fid,'%f',1)); fgetl(fid);
Conductivity = fscanf(fid,'%f',NumGridPoints); fgetl(fid);
RelativePermeability = fscanf(fid,'%f',NumGridPoints); fgetl(fid);  
Depth = fscanf(fid,'%f',NumGridPoints); fgetl(fid);

method = fscanf(fid,'%s',1); fgetl(fid);
NoiseLevel = fscanf(fid,'%f',1)/100; fgetl(fid);
FreqFromFile = round(fscanf(fid,'%f',1)); fgetl(fid); % = 1 if the frequencies are given in a file. 0 if not. 
FreqFileName = fscanf(fid,'%s',1);fgetl(fid);

if FreqFromFile
  Freq = dlmread(FreqFileName);
  Freq = Freq(:,1); % can handle when the frequencies are given in the ACPD data file.
  NumFreq = length(Freq);
else
  NumFreq = round(fscanf(fid,'%f',1)); fgetl(fid);
  FreqMin = fscanf(fid,'%f',1); FreqMax = fscanf(fid,'%f',1); 
  if NumFreq == 1
      Freq = FreqMin;
  else
      Freq = linspace(FreqMin,FreqMax,NumFreq)';  % generate a uniform set of frequency
  end
end
  
fclose(fid); % end of loading the parameters:


%---- Call the forward solver: 
Permeability = MU0*RelativePermeability; 

RefConductivity = mean(Conductivity);
RefPermeability = mean(Permeability);

% Simulate the data: 
ExSigma = PwLinCoefficient(Depth,Conductivity,RefConductivity); % Conductivity (class)
ExMu = PwLinCoefficient(Depth,Permeability,RefPermeability); % permeability (class)

PotentialDrop =  acpd1d_potentialdrop(ExSigma,ExMu,Freq,rho11,rho12,rho21,rho22,method); % measured data

% add noise to the data and smooth it out: 
PotentialDrop = real(PotentialDrop).*(1 + 2*NoiseLevel*(rand(NumFreq,1)-0.5)) + 1i*imag(PotentialDrop).*(1 + 2*NoiseLevel*(rand(NumFreq,1)-0.5));

% save the data to a file of the same format as the experimental data:
data = [Freq, real(PotentialDrop), imag(PotentialDrop)];
dlmwrite(OutputFileName,data,'delimiter',' ');
if nargin < 2
    OutputCoefFile = ['ExaCoef_',OutputFileName(1:end-4)]; % extension .mat is automatically added
end
eval(['save ' OutputCoefFile ' ExSigma ExMu data']);





