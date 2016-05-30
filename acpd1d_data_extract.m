function acpd1d_data_extract(inputfile, outputfile)
% extract only the real part and imaginary part of the potential drop data. 
% the original data file has 9 columns, we only extract three columns: 
% column 1: frequencies
% column 8: real part of the potential drop
% column 9: imaginary part of the potential drop
% inputfile: the input data file provided by John Bowler's group.
% output file: a new file with only 3 columns


data = dlmread(inputfile); % data file from Rashid
data(:,2:7) = [];

dlmwrite(outputfile,data,'delimiter',' ');
