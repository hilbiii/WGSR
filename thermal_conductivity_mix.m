function [Lamb_L,Lamb_H]=thermal_conductivity_mix(T,P,y)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [Lamb_L,Lamb_H]=thermal_conductivity_mix(T,P,y) calculates the gas
% mixture thermal conductivity at low pressure (Lamb_L) and at high 
% pressure (Lamb_H) in W/m/K. It requires the gas mixture temperature (K)
% in K, the gas mixture pressure (P) in bar and the molar fraction of all
% compoenents in the following order y=[H2,CO,CO2,H2O,N2,CH4].
% Ex: [Lamb_L,Lamb_H]=thermal_conductivity_mix(400,30,[0.1,0.1,0.1,0.6,0,0.1])
%% Parameters
%Oder of compoenents: H2, CO, CO2, H2O, N2, CH4
MW=[2.016,28.01,44.01,18.015,28.014,16.043]; % Molecular weight (kg/kmol)
Tc=[32.98,132.85,304.12,647.14,126.2,190.56]; % Critical temperature (K)
Pc=[12.93,34.94,73.74,220.64,33.98,45.99]; % Critical pressure (bar)
Vc=[64.2,93.1,94.07,55.95,90.01,98.6]; % Critical volume (cm3/mol)
Zc=[0.303,0.292,0.274,0.229,0.289,0.286]; % Critical compressibility factor(-)
DM=[0,0.1,0,1.8,0,0]; % Dipole moment (debye)
AF=[-0.217,0.045,0.225,0.344,0.037,0.011]; % Acentric factor(-)
R=8.314; % kJ/kmol/K



%% Mixture V critical
Vcc=zeros(length(Vc));
for i=1:length(Vc)
    for j=1:length(Vc)
        if i==j
            Vcc(i,j)=Vc(i);
        else
            Vcc(i,j)=(1/8)*((Vc(i)^(1/3))+(Vc(j)^(1/3)))^3;
        end
    end
end

Vcmm=zeros(length(Vc));
for i=1:length(Vc)
    for j=1:length(Vc)
        Vcmm(i,j)=y(i)*y(j)*Vcc(i,j);
    end
end
Vcm=sum(sum(Vcmm)); % cm3/mol

% Mixture T critical

Tcc=zeros(length(Tc));
for i=1:length(Tc)
    for j=1:length(Tc)
        if i==j
            Tcc(i,j)=Tc(i);
        else
            Tcc(i,j)=sqrt(Tc(i)*Tc(j));
        end
    end
end

Tcmm=zeros(length(Tc));
for i=1:length(Tc)
    for j=1:length(Tc)
        Tcmm(i,j)=y(i)*y(j)*Vcc(i,j)*Tcc(i,j);
    end
end
Tcm=sum(sum(Tcmm))/Vcm;

% Mixture AF

AFm=sum(y.*AF);

% Mixture Z critical

Zcm=0.291-0.08*AFm;

% Mixture P critical

Pcm=(Zcm*R*Tcm/Vcm)*10; % 10 is the coversion factor to transform P in bar

% Mixture MW

MWm=sum(y.*MW);


%% Thermal conductivity at low pressure

Trm=T/Tcm;

Zim=2+10.5*Trm^2;

Betham=0.7862-0.7109*AFm+1.3168*AFm^2;

[~,Cvm]=heat_capacity(T,y);

alpha=(Cvm/R)-3/2;
psi=1+alpha*((0.215+0.28288*alpha-1.061*Betham+0.26665*Zim)/(0.6366+Betham*Zim+1.061*alpha*Betham));

[Vis_L,~]=viscosity_mixture(T,P,y);

MWmn=MWm*1e-3; % conversion from kg/kmol to kg/mol
Lamb_L=3.75*psi*Vis_L*R/MWmn; % W/m/K


%% Thermal conductivity at high pressure
[Vm,~]=PRSV(T,P,y);
rho=(Vcm/Vm)*1e-3;
gama=210*(Tcm*MWm^3/(Pcm^4))^(1/6);

    if rho < 0.5
        Lamb_H=Lamb_L+(1.22e-2*(exp(0.535*rho)-1))/(gama*Zcm^5);
    elseif rho < 2
        Lamb_H=Lamb_L+(1.14e-2*(exp(0.67*rho)-1.069))/(gama*Zcm^5);
    else
        Lamb_H=Lamb_L+(2.6e-3*(exp(1.155*rho)+2.016))/(gama*Zcm^5);
    end
    
    
end