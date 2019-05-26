function j=RouletteWheelSelection(P)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
     r=rand;
    C=cumsum(P);
    j=find(r<=C,1,'first');


end

