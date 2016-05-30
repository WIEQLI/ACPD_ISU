function [Freq,PotentialDrop] = acpd1d_invprob_load_data(DataFileName)
% load the potential drop data, assign values to the vectors of frequencies
% and potential drop


dat = dlmread(DataFileName);
Freq = dat(:,1); % set of frequencies
PotentialDrop = dat(:,2) + dat(:,3)*1i; % the second column is the real part, the third column is the imaginary part of the potential drop.
