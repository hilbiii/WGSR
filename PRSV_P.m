function P=PRSV_P(T,C,y)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [Vm,Zm]=PR(T,C,y) retuns the pressure P in bar. It requires the system temperature
% (T) in K, the total concentration (C) in kmol/m3 and the molar fraction of the 
% components in the following order y=[H2,CO,CO2,H2O,N2,CH4].
% Ex: P=PRSV_P(550,0.021887597116526,[0.1,0.1,0.1,0.5,0.1,0.1])
% Algorithm based on:
% PRSV - An Improved Peng-Robinson Equation of State with New 
%  Mixing Rules for Strongly Nonideal Mixtures - Stryjek and Vera 1986
% and
% PRSV: An Improved Peng- Robinson Equation of State
% for Pure Compounds and Mixtures - Stryjek and Vera 1986

%% Parameters
% Order of the compoenents: H2, CO, CO2, H2O, N2, CH4
Tc=[32.98,132.85,304.12,647.14,126.2,190.56]; % Critical temperature (K)
Pc=[12.93,34.94,73.74,220.64,33.98,45.99]; % Critical pressure (bar)
AF=[-0.217,0.045,0.225,0.344,0.037,0.011]; % Acentric factor(-)
kappa1=[0,0,0.04285,-0.06635,0.01996,-0.00159];
R=8.314; % J/mol/K
%%
Tr=T./Tc;
kappa0=0.378893+1.4897153.*AF-0.17131848.*(AF.^2)+0.0196554.*(AF.^3);
kappa=kappa0+kappa1.*(1+Tr.^0.5).*(0.7-Tr);
alpha=(1+kappa.*(1-(Tr.^0.5))).^2;

a=0.457235*(R^2).*(Tc.^2).*alpha./(10^5.*Pc);% 10^5 is a factor to 
% convert bar to Pa. a is in J/(m^3.mol^2)
b=0.077796*R.*Tc./(10^5.*Pc); % 10^5 is a factor to 
% convert bar to Pa . b is in m3/mol

k=zeros(length(AF));
k(1,2)=0.0253;
k(1,3)=0.01202;
k(1,4)=-0.000002217*T^2-0.0009798239*T+1.9700385487;
k(1,5)=-0.036;
k(1,6)=0.2023;
k(2,1)=0.0253;
k(2,3)=-0.0314;
k(2,4)=-0.0000024279*T^2+0.0036044342*T-0.8550289761;
k(2,5)=0.0115;
k(2,6)=0.021;
k(3,1)=0.01202;
k(3,2)=-0.0314;
k(3,4)=-0.0000012681*T^2+0.0008626417*T-0.0210896419;
k(3,5)=-0.02;
k(3,6)=0.1;
k(4,1)=-0.000002217*T^2-0.0009798239*T+1.9700385487;
k(4,2)=-0.0000024279*T^2+0.0036044342*T-0.8550289761;
k(4,3)=-0.0000012681*T^2+0.0008626417*T-0.0210896419;
k(4,6)=-0.0000028361*T^2+0.0021600186*T+0.2116832245;
k(4,5)=-0.0000009013*T^2+0.0022536039*T-0.6391739866;
k(6,1)=0.2023;
k(6,2)=0.021;
k(6,3)=0.1;
k(6,4)=-0.0000028361*T^2+0.0021600186*T+0.2116832245;
k(6,5)=0.036;
k(5,1)=-0.036;
k(5,2)=0.0115;
k(5,3)=-0.02;
k(5,4)=-0.0000009013*T^2+0.0022536039*T-0.6391739866;
k(5,6)=0.036;

A=zeros(length(AF));
for i=1:length(AF)
    for j=1:length(AF)
        A(i,j)=y(i)*y(j)*(1-k(i,j))*(sqrt(a(i)*a(j)));
    end
end
Am=sum(sum(A)); % J/(m^3.mol^2)
Bm=sum(y.*b);   % m3/mol
Vm=(1/C)/1000; % m3/mol

P=(R*T/(Vm-Bm)-(Am/(Vm^2+2*Bm*Vm-Bm^2)))/10^5;

end

