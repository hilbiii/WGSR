function [Vis_Mix_Low_P,Vis_Mix_High_P]=viscosity_mixture(T,P,y)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [Vis_Mix_Low_P,Vis_Mix_High_P]=viscosity_mixture(T,P,y) requires the
% Temperature (T) in K, the pressure (P) in bar and the molar fraction (y)
% of the mixture as inputs and returns the viscosity of the mixture at low
% pressure (Vis_Mix_Low_P) and the viscosity of the mixture at high
% pressure (Vis_Mix_High_P), both in Pa.s. The molar fraction must
% follow the order y=[H2,CO,CO2,H2O,N2,CH4]
%[Vis_Mix_Low_P,Vis_Mix_High_P]=viscosity_mixture(400,30,[0.1,0.1,0.1,0.6,0,0.1])
%% Parameters
%Oder of compoenents: H2, CO, CO2, H2O, N2, CH4
MW=[2.016,28.01,44.01,18.015,28.014,16.043]; % Molecular weight (kg/kmol)
Tc=[32.98,132.85,304.12,647.14,126.2,190.56]; % Critical temperature (K)
Pc=[12.93,34.94,73.74,220.64,33.98,45.99]; % Critical pressure (bar)
Vc=[64.2,93.1,94.07,55.95,90.01,98.6]; % Critical volume (cm3/mol)
Zc=[0.303,0.292,0.274,0.229,0.289,0.286]; % Critical compressibility factor(-)
DM=[0,0.1,0,1.8,0,0]; % Dipole moment (debye)


%% Low pressure viscosity mixture calculation

Tr=T./Tc;


DDM=(52.46*DM.^2.*Pc)./Tc.^2; % Dimensionless dipole moemnt (-)

Tcm=sum(y.*Tc); % Critical temperature of the mixture (K)
Trm=T/Tcm;
Pcm=10*8.314*Tcm*(sum(y.*Zc)/sum(y.*Vc)); % Critical pressure of the mixture (bar)
MWm=sum(y.*MW); % Molecular weight of the mixture (kg/kmol)

RIV=0.176*(Tcm/((MWm^3)*(Pcm^4)))^(1/6);% Reduced inverse viscosity (1/microPoise)

Fp0=zeros(1,length(y)); % Correction factor to account for polarity
    for i = 1:length(y)
        if DDM(i) < 0.0222
            Fp0(i)=1;
        elseif DDM(i) < 0.075
            Fp0(i)=1+30.55*(0.292-Zc(i))^(1.72);
        else
            Fp0(i)=1+30.55*((0.292-Zc(i))^(1.72))*abs(0.96+0.1*(Tr(i)-0.7));
        end
    end
    
Fq0=ones(1,length(y));% Correction factor to account for quantum effects. Only for H2, He and D2

    if (Tr(1)-12) < 0
        fac=-1;
    else
        fac=1;
    end
    
Fq0(1)=1.22*(0.76)^0.15*(1+0.00385*(((Tr(1)-12)^2)^(1/MW(1)))*fac);

% The A factor bellow is calculated using the value of MW, but sometimes 
% one of the componets might not be present in the mixture. Therefore, only
% the MW of the present compoents must be taken into account.
MWa=nonzeros(y.*MW);
yc=nonzeros(y);
MWc=MWa./yc;
    if (max(MWc)/min(MWc)) > 9
        A=1-0.01*(max(MWc)/min(MWc))^0.87; % Factor for quantum effects corrections
    else
        A=1;
    end
%%    
Fp0m=sum(y.*Fp0);
Fq0m=sum(y.*Fq0)*A;

Z1=(0.807*(Trm^0.618)-0.357*exp(-0.449*Trm)+0.34*exp(-4.058*Trm)+0.018)*Fq0m*Fp0m;

Vis_Mix_Low_P=(Z1/RIV)*1e-7; % (Pa.s)

%% High pressure viscosity mixture calculation

Prm=P/Pcm; % Reduced pressure of the mixture
    if Trm < 1
        alpha_vis=3.262+14.98*(Prm)^5.508;
        betha_vis=1.39+5.746*Prm;
        Z2=0.6+0.76*Prm^alpha_vis+((6.990*Prm^betha_vis)-0.6)*(1-Trm);
    else
        a_vis=((1.245e-3)/Trm)*exp(5.1726*Trm^-0.3286);
        b_vis=a_vis*(1.6553*Trm-1.2723);
        c_vis=(0.4489/Trm)*exp(3.0578*Trm^-37.7332);
        d_vis=(1.7368/Trm)*exp(2.231*Trm^-7.6351);
        e_vis=1.3088;
        f_vis=0.9425*exp(-0.1853*Trm^0.4489);
        Z2=Z1*(1+((a_vis*Prm^e_vis)/((b_vis*Prm^f_vis)+(1+c_vis*Prm^d_vis)^-1)));
    end
    
Y=Z2/Z1;
Fpm=(1+(Fp0m-1)*Y^-3)/Fp0m;
Fqm=(1+(Fq0m-1)*((Y^-1)-0.007*(log(Y))^4))/Fq0m;

Vis_Mix_High_P=(Z2*Fpm*Fqm/RIV)*1e-7; % (Pa.s)
end