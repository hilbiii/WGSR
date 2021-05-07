function [D_p,D_f]=diffusion_c(T,P,y,C)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function calculates the effective diffusion coeeficient of a gas mixture
% T must be in K, P in bar and y is the molar fraction vector 
% (y=[H2,CO,CO2,H2O,N2,CH4]) D_p is the intraparticle diffusion (m^2/s)
% D_f is the bulk fluid diffusion (m^2/s)
% Ex: [D_p D_f]=diffusion_c(500,10,[0.1 0.1 0.2 0.2 0.1 0.3],1)
% C is a variable used to consider or not Knudsen diffusion. If C=1, 
% Knudesen diffusion is considered. Otherwise, it's not. Default C=1.

%% Parameters
dpore=9e-7;%cm
% Order of compoenents: H2, CO, CO2, H2O, N2, CH4
MW=[2.016,28.01,44.01,18.015,28.014,16.043]; % Molecular weight (kg/kmol)
DV=[6.12 18.9 26.9 12.7 17.9 24.14];% Diffusion volumes 

%% 
% Default C=1.
if nargin < 4
C=1;
end


%%
% Diffusion of the ith coponent into the jth component
Dij=zeros(length(y),length(y)); 
% Molar fraction of the jth component on a ith component free basis
yf=zeros(length(y),length(y)); 

Dimc=zeros(length(y),length(y)); 
for i=1:length(y)
    for j=1:length(y)
        if i~=j
            if y(i)==0
                Dij(i,j)=0;
                yf(i,j)=0;
                Dimc(i,j)=0;
            else
                Dij(i,j)=(1e-3*(T^1.75)*sqrt((1/MW(i))+(1/MW(j))))/...
                    ((P/1.01325)*(((DV(i)^(1/3))+(DV(j)^(1/3)))^2));
                yf(i,j)=y(j)/(1-y(i));
                Dimc(i,j)=(yf(i,j)/Dij(i,j));
            end
        end
     end
end
% Diffusion coefficient of the ith component into the gas mixture (cm2/s)
Dim=1./sum(Dimc');

%%
% Knudsen diffusion
Dk=zeros(1,length(y)); 
for i=1:length(y)
    if y(i)==0
        Dim(i)=0;
        Dk(i)=0;
    else
        Dim(i)=Dim(i);
        Dk(i)=4850*dpore*sqrt(T/MW(i));
    end
end

%%
% Effective diffusion coefficient of the ith component (m2/s)
D_p=zeros(1,length(y));
if C==1
    for i=1:length(y)
        D_p(i)=1e-4*0.11*(1/((1/Dk(i))+(1/Dim(i))));
    end
else
    for i=1:length(y)
        D_p(i)=1e-4*0.11*Dim(i);
    end
end
%%
% Diffusion coefficient of the ith component into the gas mixture inside
% the reactor (m2/s).
D_f=1e-4*Dim;
end