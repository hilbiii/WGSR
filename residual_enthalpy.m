function H=residual_enthalpy(T,P,y)
%It returns the Residual Enthalpy (J/mol)
% Fluid total V(m3/kmol)
[Vm,~]=PRSV(T,P,y);
R=8.314;
% Parameters
% Order of the compoenents: H2, CO, CO2, H2O, N2, CH4
Tc=[32.98,132.85,304.12,647.14,126.2,190.56]; % Critical temperature (K)
Pc=[12.93,34.94,73.74,220.64,33.98,45.99]; % Critical pressure (bar)
AF=[-0.217,0.045,0.225,0.344,0.037,0.011]; % Acentric factor(-)
kappa1=[0,0,0.04285,-0.06635,0.01996,-0.00159];

% Binary interaction parameters kij
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

% Derivative of the binary interaction parameter dkij
dk=zeros(length(AF));

dk(1,4)=-2*0.000002217*T-0.0009798239;
dk(2,4)=-2*0.0000024279*T+0.0036044342;
dk(3,4)=-2*0.0000012681*T+0.0008626417;
dk(4,1)=-2*0.000002217*T-0.0009798239;
dk(4,2)=-2*0.0000024279*T+0.0036044342;
dk(4,3)=-2*0.0000012681*T+0.0008626417;
dk(4,5)=-2*0.0000009013*T+0.0022536039;
dk(4,6)=-2*0.0000028361*T+0.0021600186;
dk(5,4)=-2*0.0000009013*T+0.0022536039;
dk(6,4)=-2*0.0000028361*T+0.0021600186;

% PRSV parameters
Tr=T./Tc; % Reduced temperature
kappa0=0.378893+1.4897153.*AF-0.17131848.*(AF.^2)+0.0196554.*(AF.^3);
kappa=kappa0+kappa1.*(1+Tr.^0.5).*(0.7-Tr);
alpha=(1+kappa.*(1-(Tr.^0.5))).^2;
a=0.457235*(R^2).*(Tc.^2).*alpha./(Pc.*10^2); % Pc*10^2 is a coversion 
% factor from bar to kJ/m3
bi=0.077796*R.*Tc./(Pc.*10^2); 
b=sum(y.*bi);

A=zeros(length(AF));
for i=1:length(AF)
    for j=1:length(AF)
        A(i,j)=y(i)*y(j)*(1-k(i,j))*(sqrt(a(i)*a(j)));
    end
end

X=sum(sum(A)); % "X" is "a" in Eq.1 of the papper " PRSV: An Improved 
% Peng- Robinson Equation of State for Pure Compounds
% and Mixtures - Stryjek and Vera 1986"

% Derivative dadT of the ith component
da=2*(0.457235*(8.314^2)*(Tc.^2)./(Pc*10^2)).*(1+kappa0+0.7*kappa1-...
    1.7*kappa1.*Tr-kappa0.*(Tr.^0.5)+kappa1.*(Tr.^2)).*...
    (-1.7*kappa1./Tc-0.5*kappa0./((T*Tc).^0.5)+2*kappa1*T./Tc.^2);

% Derivative of a for the mixture
dA=zeros(length(AF));
for i=1:length(AF)
    for j=1:length(AF)
        if i~=4 && j~=4
            dA(i,j)=(y(i)*y(j)*0.5*(1-k(i,j))/(sqrt(a(i)*a(j))))*...
            (da(i)*a(j)+a(i)*da(j));
        else
            dA(i,j)=y(i)*y(j)*(((0.5/(sqrt(a(i)*a(j))))*...
                (da(i)*a(j)+a(i)*da(j)))-(dk(i,j)*sqrt(a(i)*a(j))+...
                (0.5*k(i,j)/sqrt(a(i)*a(j)))*(da(i)*a(j)+a(i)*da(j))));
        end
    end
end
dX=sum(sum(dA));

H=((X-T*dX)/(2*sqrt(2)*b))*log((Vm+b*(1-sqrt(2)))/(Vm+b*(1+sqrt(2))))+...
    P*100*Vm-R*T;

end