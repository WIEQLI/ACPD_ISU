function [rho11,rho12,rho21,rho22,PerPar,Nz,Depth,RegPar] = acpd1d_invprob_lin_get_parameters(ParameterFile)
% get parameters from the parameter file for the linearized ACPD inverse problem. 
% Input: ParameterFile: a structured text file containing the following
% rows:
%---
% % 0.0015  0.003 0.003 0.0015          Distance_between_pins_rho11_rho12_rho21_rho22
% % 0.3                                 Perturbation parameter "epsilon" in the linearized model
% % 15                                  Number_of_grid_points_in_depth:                   
% % -0.002 0                            interval of depth where the parameters are estimated. the second parameter must always be zero! 
% % 1                                   Regularization parameters
% ---
% OUTPUT: 
% rho11, rho12, rho21, rho22: distances between the pins
% PerPar: perturbation parameter, epsilon, in the linearized model
% Nz: number of grid points in depth.
% Depth: a vector of grid points in depth. 
% RegPar: regularization parameter for the linearized problem. 
% 
% Updated on March 25, 2016.


% ---- open the input file: 
fid = fopen(ParameterFile);
if fid == -1
    error('The parameter file does not exist, stopped.');
end

rho11 = fscanf(fid,'%f',1); rho12 = fscanf(fid,'%f',1);
rho21 = fscanf(fid,'%f',1); rho22 = fscanf(fid,'%f',1); 
fgetl(fid); % get the rest of the first row.

PerPar = fscanf(fid,'%f',1); fgetl(fid); % perturbation parameter 'epsilon'
Nz = round(fscanf(fid,'%f',1)); fgetl(fid); % number of points in depth where the coefficients to be estimated.

if Nz < 2
  error('The number of points in depth must be at least 2. Please modify the input file');
end

Zmin = fscanf(fid,'%f',1); Zmax = fscanf(fid,'%f',1); fgetl(fid);
Depth = linspace(Zmin,Zmax,Nz); % a uniform partition in depth.
RegPar = fscanf(fid,'%f',1);  % regularization parameter.

fclose(fid); % end of loading the parameters:
