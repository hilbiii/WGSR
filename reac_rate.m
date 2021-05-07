function r=reac_rate(Tf,Pf,yf,v,dp)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function calculates the actual reaction rate, considering mass transfer
% and diffusion inside the catalytic particle. T must be in K, P in bar, 
% y is the molar fraction vector [H2, CO, CO2, H2O, N2, CH4] 
% v is the fluid velocity (m/s) and dp is the particle diameter (m)
% PRSV e C calculado pelo somatório de Ci
% Ex.: r=reac_rate(450,1.01325,[0.35 0.1 0.11 0.4 0.03 0.01],1,2e-3)
%% 
% Transfering input variables
p.Pf=Pf;            % bar
p.Tf=Tf;            % K
p.yf=yf;            % -
p.Rg=8.314;         % J/mol/K
p.R=dp/2;           % m  
p.lamb=0.43;        % W/m/K
p.v=v;              % m/s

% Setting reference values
[Vm,~]=PRSV(p.Tf,p.Pf,p.yf);
p.C_ref=1/Vm;   % kmol/m3
p.T_ref=500;    % K

% Setting orthogonal collocation matrices
p.N=10; % Number of internal collocation points
[x,A,B] = OCMatrices(p.N,0,0);

% Setting intial guesses
X0=0.1*ones(7,p.N+2);   % y
X0(7,:)=10*X0(7,:);     % T

% Setting solver options
op=optimoptions('fsolve','TolFun',1e-15,...
'TolX',1e-15,'MaxFunEvals',10e6,...
'MaxIter',1e6,'Display','off');

X = fsolve(@Balances,X0,op,p,x,A,B);

% Mass (m/s) transfer coefficient
[K_m,~]=massheatpartcoef(p.Tf,p.Pf,p.yf,p.v,dp);

% Fluid total concentration Cf at bulk fluid condition
[Vm,~]=PRSV(p.Tf,p.Pf,p.yf);
Cf=1/Vm;
%% Real reaction rate
% CO is consumed. Therefore its rate is negative. To change that, a minus
% sign is placed in front of the equation. This equation means that the
% amount of CO transpoted through the risitance film that involves the
% catalyst particle must be equal to the rate of consuption of CO.

r=-(3/p.R)*(K_m(2))*(p.C_ref*X(2,end)-Cf*p.yf(2))*1000; %(mol/m^3/s)

%% Figures
rp=x*p.R*100;
subplot(1,2,1);
plot(rp,p.C_ref.*X(1,:),'k-',rp,p.C_ref.*X(2,:),'k.',rp,p.C_ref.*X(3,:)...
    ,'k.-',rp,p.C_ref*X(4,:),'k-.',rp,p.C_ref*X(5,:),'k.-.',...
    rp,p.C_ref*X(6,:),'k+',rp(end),p.C_ref*p.yf(1),'ko',...
    rp(end),p.C_ref*p.yf(2),'ko',rp(end),p.C_ref*p.yf(3),'ko'...
    ,rp(end),p.C_ref*p.yf(4),'ko',rp(end),p.C_ref*p.yf(5),'ko',...
    rp(end),p.C_ref*p.yf(6),'ko');
legend('H2','CO','CO_2','H_2O','N_2','CH_4','Location','best')
xlabel('Particle radius (cm)')
ylabel('Concentration (kmol/m^3)')
xlim([0 (rp(end)+0.01)])
ylim([0 0.013])
txt1 = '(a)';
text(0.01,0.002,txt1)

subplot(1,2,2)
plot(rp,p.T_ref*X(7,:),'k-',rp(end),p.Tf,'ko')
xlabel('Particle radius (cm)')
ylabel('Temperature (K)')
xlim([0 (rp(end)+0.01)])
txt2 = '(b)';
text(0.01,450.03,txt2)
end

function F=Balances(X0,p,x,A,B)
% global C P y
% Transfering the initial guesses to the respective variables
C_H2=X0(1,:);
C_CO=X0(2,:);
C_CO2=X0(3,:);
C_H2O=X0(4,:);
C_N2=X0(5,:);
C_CH4=X0(6,:);
T=X0(7,:);

%% At x=0. 
% Total dimentionless fluid concentration
C=zeros(1,p.N+2);
C(1)=C_H2(1)+C_CO(1)+C_CO2(1)+C_H2O(1)+C_N2(1)+C_CH4(1);

% Molar fraction
y=zeros(p.N+2,length(p.yf));
y(1,1)=C_H2(1)/C(1);    % H2
y(1,2)=C_CO(1)/C(1);    % CO
y(1,3)=C_CO2(1)/C(1);   % CO2
y(1,4)=C_H2O(1)/C(1);   % H2O
y(1,5)=C_N2(1)/C(1);    % N2
y(1,6)=C_CH4(1)/C(1);   % CH4

% Boundary condition
F(1,1)=A(1,:)*C_H2';
F(2,1)=A(1,:)*C_CO';
F(3,1)=A(1,:)*C_CO2';
F(4,1)=A(1,:)*C_H2O';
F(5,1)=A(1,:)*C_N2';
F(6,1)=A(1,:)*C_CH4';
F(7,1)=A(1,:)*T';

%% 0<x<1
for i=2:p.N+1
    
    % Total dimentionless fluid concentration
    C(i)=C_H2(i)+C_CO(i)+C_CO2(i)+C_H2O(i)+C_N2(i)+C_CH4(i);

    % Molar fraction
    y(i,1)=C_H2(i)/C(i);    % H2
    y(i,2)=C_CO(i)/C(i);    % CO
    y(i,3)=C_CO2(i)/C(i);   % CO2
    y(i,4)=C_H2O(i)/C(i);   % H2O
    y(i,5)=C_N2(i)/C(i);    % N2
    y(i,6)=C_CH4(i)/C(i);   % CH4
    
    % Rate constant (kmol/m3/s)
    k=5904*(1-0.55)*2.96e5*exp(-47400/(p.Rg*p.T_ref*T(i)));
    
    % Equilibrium constant
    K=1/(exp((13.148)-(5639.5/(p.T_ref*T(i)))-1.077*log(p.T_ref*T(i))-...
        (5.44*10^-4)*p.T_ref*T(i)-(1.125*10^-7)*((p.T_ref*T(i))^2)+...
        (49170/((p.T_ref*T(i))^2))));
    
    % Heat of reaction (J/kmol)
    mDH=-1e6*(-47.617+1.302e-2*T(i)*p.T_ref-0.126e-5*((T(i)*p.T_ref)^2)+...
        0.791e3*((T(i)*p.T_ref)^-1));
    
    % Raction rate (kmol/m3/s)
    r=((k*(y(i,2)*y(i,4)-y(i,3)*y(i,1)/K)))/3600;
    
    % Diffusion coeefficient of the ith component on the gas mixture Dim 
    % (m^2/s)
    [D,~]=diffusion_c(T(i)*p.T_ref,p.Pf,y(i,:),1);
    
    %Molar Balance 
    % H2 molar balance    
    F(1,i)=(2/x(i))*(A(i,:)*C_H2')+B(i,:)*C_H2'+(p.R^2)*r/(D(1)*p.C_ref);
    
    % CO molar balance
    F(2,i)=(2/x(i))*(A(i,:)*C_CO')+B(i,:)*C_CO'-(p.R^2)*r/(D(2)*p.C_ref);
    
    % CO2 molar balance
    F(3,i)=(2/x(i))*(A(i,:)*C_CO2')+B(i,:)*C_CO2'+(p.R^2)*r/(D(3)*p.C_ref);
    
    % H2O molar balance
    F(4,i)=(2/x(i))*(A(i,:)*C_H2O')+B(i,:)*C_H2O'-(p.R^2)*r/(D(4)*p.C_ref);
    
    % N2 molar balance
    F(5,i)=(2/x(i))*(A(i,:)*C_N2')+B(i,:)*C_N2';
    
    % CH4 molar balance
    F(6,i)=(2/x(i))*(A(i,:)*C_CH4')+B(i,:)*C_CH4';
    
    % Energy balance    
    F(7,i)=(2/x(i))*A(i,:)*T'+B(i,:)*T'+mDH*r*p.R^2/(p.lamb*p.T_ref);
    
end
%% At x=1
% Total dimentionless fluid concentration
C(p.N+2)=C_H2(p.N+2)+C_CO(p.N+2)+C_CO2(p.N+2)+...
    C_H2O(p.N+2)+C_N2(p.N+2)+C_CH4(p.N+2);

% Molar fraction
y(p.N+2,1)=C_H2(p.N+2)/C(p.N+2);    % H2
y(p.N+2,2)=C_CO(p.N+2)/C(p.N+2);    % CO
y(p.N+2,3)=C_CO2(p.N+2)/C(p.N+2);   % CO2
y(p.N+2,4)=C_H2O(p.N+2)/C(p.N+2);   % H2O
y(p.N+2,5)=C_N2(p.N+2)/C(p.N+2);    % N2
y(p.N+2,6)=C_CH4(p.N+2)/C(p.N+2);   % CH4

% Mass (m/s) and heat (W/m^2/K) transfer coefficients
[K_m,h]=massheatpartcoef(p.Tf,p.Pf,p.yf,p.v,2*p.R);

% Diffusion coeefficient of the ith compoent on the gas mixture Dim 
% (m^2/s)
[D,~]=diffusion_c(T(p.N+2)*p.T_ref,p.Pf,y(p.N+2,:),1);

% Fluid total concentration Cf at bulk fluid condition
[Vm,~]=PRSV(p.Tf,p.Pf,p.yf);
Cf=1/Vm;

% Boundary condition
% H2
F(1,p.N+2)=A(p.N+2,:)*C_H2'+...
        (K_m(1)*p.R/D(1))*(C_H2(p.N+2)-Cf*p.yf(1)/p.C_ref);

% CO
F(2,p.N+2)=A(p.N+2,:)*C_CO'+...
        (K_m(2)*p.R/D(2))*(C_CO(p.N+2)-Cf*p.yf(2)/p.C_ref);

% CO2
F(3,p.N+2)=A(p.N+2,:)*C_CO2'+...
        (K_m(3)*p.R/D(3))*(C_CO2(p.N+2)-Cf*p.yf(3)/p.C_ref);

% H2O
F(4,p.N+2)=A(p.N+2,:)*C_H2O'+...
        (K_m(4)*p.R/D(4))*(C_H2O(p.N+2)-Cf*p.yf(4)/p.C_ref);
    
% N2    
F(5,p.N+2)=A(p.N+2,:)*C_N2'+...
        (K_m(5)*p.R/D(5))*(C_N2(p.N+2)-Cf*p.yf(5)/p.C_ref);
% CH4    
F(6,p.N+2)=A(p.N+2,:)*C_CH4'+...
        (K_m(6)*p.R/D(6))*(C_CH4(p.N+2)-Cf*p.yf(6)/p.C_ref);
    
% T
F(7,p.N+2)=A(p.N+2,:)*T'+(h*p.R/p.lamb)*(T(p.N+2)-p.Tf/p.T_ref);

end