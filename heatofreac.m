function mDH=heatofreac(T,P,y)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DHm=heatofreac(T,P,y) requires the temperature (T) in K, the 
% pressure (P) in bard and the components molar fractions (y)
% It returns heat of reaction (kJ/kmol) 
% The reference temperature is 298.15 K. 
% y=[H2 CO CO2 H2O N2 CH4]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Standard Enthalpy of Formation
R=8.314; % Ideal gas constant (J/mol/K)
T0=298.15; % Reference temperature (K)
h0=[0,-110530,-393510,-241810]; % Standard enthalpy of formation (kJ/kmol)
H0=h0(1)+h0(3)-h0(2)-h0(4); % Heat of reaction at T0 and 1 bar

%% Ideal Gas Enthalpy
a0=[2.883 3.912 3.259 4.395];
a1=[3.681 -3.913 1.356 -4.186]*1e-3;
a2=[-0.772 1.182 1.502 1.405]*1e-5;
a3=[0.692 -1.302 -2.374 -1.564]*1e-8;
a4=[-0.213 0.515 1.056 0.632]*1e-11;

Da0=a0(1)+a0(3)-a0(2)-a0(4);
Da1=a1(1)+a1(3)-a1(2)-a1(4);
Da2=a2(1)+a2(3)-a2(2)-a2(4);
Da3=a3(1)+a3(3)-a3(2)-a3(4);
Da4=a4(1)+a4(3)-a4(2)-a4(4);
Hig=R*(Da0*(T-T0)+(Da1/2)*(T^2-T0^2)+(Da2/3)*(T^3-T0^3)+...
    (Da3/4)*(T^4-T0^4)+(Da4/4)*(T^5-T0^5)); % kJ/kmol
y=y';
%% Residual Enthalpy 
Hres=residual_enthalpy(T,P,y);

%% Enthalpy of Reaction
mDH=-(H0+Hig+Hres);

end


