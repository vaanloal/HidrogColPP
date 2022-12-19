function [mu_xil] = mu_xil(T)
T = T + 273.15;

A=	-4.617;
B=	775.41;
C=	-0.94611;
D=	0;
E=	0;
F=	0;
G=	0;



mu_xil = exp(A + B./T + C* log(T)) * 1000; % cP
end