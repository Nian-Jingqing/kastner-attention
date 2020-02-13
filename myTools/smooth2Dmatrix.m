function out = smooth2Dmatrix(in, sigma)
% this function takes in a 2D matrix, performs smoothing and returns the output as a matrix

in_struct(1).field = in;
r_struct = R.Rstruct(in_struct);
r_struct.smoothFieldInR('field', 'smoothed', sigma, 1);
out = r_struct.r(1).smoothed;
%