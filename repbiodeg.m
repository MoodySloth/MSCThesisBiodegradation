function [R,B,M,S] = repbiodeg(Blank,Mat,T);
R = (Mat - Blank);
B = (R./T)*100;
M = nanmean(B,2);
S = nanstd(B,0,2);
end