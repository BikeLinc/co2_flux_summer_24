function [mol] = PPM_TO_MOL(ppm)
%PPM_TO_MOL Summary of this function goes here
%   Detailed explanation goes here
    mol = (ppm*P)/(1e6*8.3145*T);
end

