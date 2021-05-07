function [Cpm,Cvm]=heat_capacity(T,y)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [Cpm,Cvm]=heat_capacity(T,y) requires the temperature (T) in K and the
% component molar fraction in the following order y=[H2,CO,CO2,H2O,N2,CH4].
% It returns the mixture heat capacity at constant pressure (Cpm) and at
% constant volume (Cvm). Both in kJ/kmol/K. The reference temperature is 298.15 K. 
%ex.: [Cpm,Cvm]=heat_capacity(550,[0.1,0.1,0.1,0.5,0.1,0.1])
R=8.314; %kJ/kmol/K
a0=[2.883 3.912 3.259 4.395 3.539 4.568];
a1=[3.681 -3.913 1.356 -4.186 -0.261 -8.975]*1e-3;
a2=[-0.772 1.182 1.502 1.405 0.007 3.631]*1e-5;
a3=[0.692 -1.302 -2.374 -1.564 0.157 -3.407]*1e-8;
a4=[-0.213 0.515 1.056 0.632 -0.099 1.091]*1e-11;

Cpm=sum(R*y.*(a0+a1*T+a2*T^2+a3*T^3+a4*T^4));
Cvm=Cpm-R;

end