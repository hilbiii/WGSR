function K=equiconst(T,P,y)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% K=equiconst(T,P,y) requires the temperature (T) in K, the 
% pressure (P) in bard and the components molar fractions (y)
% It returns the equilibrium constant 
% The reference temperature is 298.15 K. 
% y=[H2 CO CO2 H2O N2 CH4]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lnK0=11.804666459469091; % ln of K at 298.15 K and 1 bar
T0=298.15;
op = odeset('RelTol',1e-4,'AbsTol',1e-4);
R=8.314; % kJ/kmol/K
[T,lnK] = ode45(@derivate,[T0 T],lnK0,op,P,y,R);

K=exp(lnK(end));

end

function dlnK=derivate(T,K,P,y,R)

DH=-(heatofreac(T,P,y));

dlnK=DH/(R*T^2);


end