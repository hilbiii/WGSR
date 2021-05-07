function [K_m,h]=massheatpartcoef(T,P,y,v,dp)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [K_m,h]=massheatpartcoef(T,P,y,v,dp) requires the
% Temperature (T) in K, the pressure (P) in bar, the molar fraction (y)
% of the mixture, the fuid velocity (m/s) and the particle diameter (m) 
% as inputs and returns the mass (K_m)(m/s) transefr coefficients of the
% ith component and heat (h)(W/m2.K) transfer coefficient.
% The molar fraction must follow the order y=[H2,CO,CO2,H2O,N2,CH4]
% Ex: [K_m,h]=massheatpartcoef(500,2,[0.2,0.1,0.2,0.3,0.1,0.1],1,0.01)

%%
MW=[2.016,28.01,44.01,18.015,28.014,16.043]; % Molecular weight (kg/kmol)
[~,Diff]=diffusion_c(T,P,y);
[V,~]=PRSV(T,P,y);
[~,Vis_Mix]=viscosity_mixture(T,P,y);
Cpm=cp_real(T,P,y);
[~,Lamb]=thermal_conductivity_mix(T,P,y);
rho=sum(y.*MW)/V; % Mixture density (kg/m3)

%% Mass transfer coefficient
Re=rho*dp*v/Vis_Mix; % Reynolds number
Sc=zeros(1,length(Diff));
for i=1:length(Diff)
    if Diff(i)==0
        Sc(i)=0; % Schmidt number
    else
        Sc(i)=Vis_Mix/(rho*Diff(i)); % Schmidt number
    end
end
Sh=2+1.1*Re^(0.6).*Sc.^(0.33);% Sherwood number
K_m=Sh.*Diff/dp; % Mass tranfer coefficient (m/s)

%% Heat transfer coefficient between particle and bulk fluid
Pr=Cpm*Vis_Mix/(Lamb); % Prandtl number
Nu=2+1.1*Re^(0.6)*Pr^(0.33); % Nusselt number
h=Nu*Lamb/dp; % Heat transfer coeficient (W/m2.K)

end