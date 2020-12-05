function [R,BSm,M,S] = repbiodegsm(Blank,Mat,T);
% Blank = Matrix of blank vessel values
% Mat = Matrix of the accumulated values of a material
% T = Theoretical CO2 of Mat
windowSize = 10;
b = (1/windowSize)*ones(1,windowSize);
a = 1;

R = (Mat - Blank);
B = (R./T)*100;
BSm = filter(b,a,B);
M = nanmean(BSm,2);
S = nanstd(BSm,0,2);
end

