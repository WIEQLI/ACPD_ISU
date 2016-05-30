function output = integral(f,z1,z2,tolerance)
% make this function consistent with the newer version

output = quad(f,z1,z2,tolerance);
