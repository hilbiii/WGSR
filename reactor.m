function [P_out,y_out,L_inert,L]=reactor(T,P,n,y,dp,dt,T_inert)
% Reactor tube model. Ir returns the reactor outlet pressure, composition,
% length and inert layer length if any.,

global check count z_in z_out 
count=0;
z_in=[];
z_out=[];
check=0;
%%
% Effectiviness factor equation coefficients.
load('coef_dp1.mat')
load('coef_dp2.mat')
load('coef_dp3.mat')
load('coef_dp4.mat')
load('coef_dp5.mat')
p.dp1=coef_dp1;
p.dp2=coef_dp2;
p.dp3=coef_dp3;
p.dp4=coef_dp4;
p.dp5=coef_dp5;

p.P_in=P;     % Presure at reactor inlet (bar)
p.T_in=T;      % Temperture at reactor inlet (K)
p.y_in=y; % Component molar fraction at reactor inlet
p.T_inert_cool=T_inert; % Temperature to be achived in the inert layer (K)
p.dt=dt;   % Reactor tube diameter (m) 
p.Rr=p.dt/2;       % Reactor tube radius (m) 
p.At=pi()*p.Rr^2;   % Reactor tube cross sectional area (m2)
p.Rg=8.314;         % Ideal gas constant (J/mol/K)
p.L=10;         % Reactor tube length (m)
p.n_in=n;      %  mol/s

p.Tcf=T_inert-5;          % Cooling fluid temperature (K)
p.U=35;             % Overall heat transfer coefficient (W/m2/K)

[Vmin,~]=PRSV(p.T_in,p.P_in,p.y_in);  % PRSV
p.C_in=(1/Vmin)*1000;     % Total fluid concentration at inlet (mol/m3)
p.Q=p.n_in/p.C_in;            % Fluid volumetric flow rate at inlet (m3/s)
p.v_in=p.Q/p.At;               % m/s
p.MW=[2.016,28.01,44.01,18.015,28.014,16.043]./1000;       % Molecular
                                                           % weight(kg/mol)
p.N_var=7;      % Number of variables
p.NC=6;         % Number of components
p.dp=dp;   % Particle diamter (m)
p.lamb_cat=0.43;    % Catalyst particle thermal conductivity (W/m/K)

%% Setting orthogonal collocation matrices
p.NR=18; % Number of internal collocation points in the radial coordinate
[xr,AR,BR] = OCMatrices(p.NR,0,0);

%% Setting reference variables
p.v_ref=0.2;                 % m/s
p.C_ref=0.5*p.C_in;          % mol/m3
p.T_ref=500;                 % K
p.D_ref=1e-5;                % m2/s

%% Creating variables vectors
% Rate constant (kmol/m3/s)
p.kr=zeros((p.NR+2),1);

% Equilibrium constant
p.K=zeros((p.NR+2),1);

% Fluid total concentration (dimentionless)
p.Ct=zeros((p.NR+2),1);

% Molar Fraction
p.y=zeros(p.NC,(p.NR+2));

% Diffusion coefficient of th ith component on the mixture
p.Dmi=zeros(p.NC,(p.NR+2));

% Diffucion coefficient axial coordinate
p.Dr=zeros(p.NC,(p.NR+2));

% Fluid viscosity (Pa.s)
p.mu=zeros((p.NR+2),1);

% Fluid thermal conductivity (W/m/K)
p.lambF=zeros((p.NR+2),1);

% Fluid specific heat (J/mol/K)
p.Cp=zeros((p.NR+2),1);

% Reynolds number
p.Re=zeros((p.NR+2),1);

% Prandtl number
p.Pr=zeros((p.NR+2),1);

% Thermal conductivity - Static contribuition
p.lamb0=zeros((p.NR+2),1);

% Thermal conductivity - Dinamic contribuition radial
p.lambDr=zeros((p.NR+2),1);

% Thermal conductivity - Effective radial
p.lambR=zeros((p.NR+2),1);

% dCdT at constant pressure
p.dC=zeros((p.NR+2),1);

% Heat of reaction (J/mol)
p.mDH=zeros((p.NR+2),1);

% Raction rate (mol/m3/s)
p.r=zeros((p.NR+2),1);

%% Calculating the velocity profile and pressure drop
[p.v,p.dPdz]=momentum(p.T_in,p.P_in,p.y_in,p.Q,p.dp,p.dt,p.NR);

p.P_drop=p.L*p.dPdz/1e5; % Total pressure drop (bar)
                                                   
% Porosity as a function of the radial coodinate
p.epsilon=zeros(1,p.NR+2);
for i=1:p.NR+2
    p.eps(i)=0.4*(1+1.36*exp(-5*(p.Rr-xr(i)*p.Rr)/p.dp));
end

% Setting Mass matrix. It specificies which equations are algebraic (0) and
% which ones are differential (1). Only the main diagonal is used.
M_init =zeros(p.N_var*(p.NR+2)+1,p.N_var*(p.NR+2)+1);
for k=1:p.N_var
    for i=(k-1)*(p.NR+2)+1:k*(p.NR+2)
        for j=(k-1)*(p.NR+2)+1:k*(p.NR+2)
            if i==j
               M_init(i,j)=1;
            end
        end
    end
end
vec_init=zeros((p.NR+2),1);
for i=2:(p.NR+1)
       vec_init(i,1)=1;
end
for k=1:p.N_var
   vec((k-1)*(p.NR+2)+1:(k)*(p.NR+2),1)=vec_init(:,1);
end
M=zeros(p.N_var*(p.NR+2)+1,p.N_var*(p.NR+2)+1);
for i=1:p.N_var*(p.NR+2)
   M(1:(p.N_var*(p.NR+2)),i)=M_init(i,1:(p.N_var*(p.NR+2)))'.*vec;
end
M(end,end)=1;

% Setting the velocity to be dimentionless
p.v=(p.v./p.v_ref);

% Integration length (dimentionless)
z=[0 1];

% Initial condition
x0=ones(1,p.N_var*(p.NR+2)+1);
x0(1:(p.NR+2))=x0(1:(p.NR+2)).*p.C_in*p.y_in(1)./p.C_ref;
x0((p.NR+2)+1:2*(p.NR+2))=x0((p.NR+2)+1:2*(p.NR+2)).*p.C_in...
    *p.y_in(2)./p.C_ref;
x0(2*(p.NR+2)+1:3*(p.NR+2))=x0(2*(p.NR+2)+1:3*(p.NR+2)).*p.C_in...
    *p.y_in(3)./p.C_ref;
x0(3*(p.NR+2)+1:4*(p.NR+2))=x0(3*(p.NR+2)+1:4*(p.NR+2)).*p.C_in...
    *p.y_in(4)./p.C_ref;
x0(4*(p.NR+2)+1:5*(p.NR+2))=x0(4*(p.NR+2)+1:5*(p.NR+2)).*p.C_in...
    *p.y_in(5)./p.C_ref;
x0(5*(p.NR+2)+1:6*(p.NR+2))=x0(5*(p.NR+2)+1:6*(p.NR+2))*p.C_in...
    *p.y_in(6)./p.C_ref;
x0(6*(p.NR+2)+1:7*(p.NR+2))=x0(6*(p.NR+2)+1:7*(p.NR+2))*p.T_in./p.T_ref;
x0(1,end)=0; % Heat Flux Initial Condition


% Setting solver options
op = odeset('Mass',M,'RelTol',1e-6,'AbsTol',1e-8,'Events',@events);

% Calling solver for integration
[z,C] = ode15s(@diff,z,x0,op,p,xr,AR,BR);


for i=1:length(z)
    for j=1:p.NR+2
        C_H2(j,i)=C(i,0*(p.NR+2)+j)*p.C_ref;
        C_CO(j,i)=C(i,1*(p.NR+2)+j)*p.C_ref;
        C_CO2(j,i)=C(i,2*(p.NR+2)+j)*p.C_ref;
        C_H2O(j,i)=C(i,3*(p.NR+2)+j)*p.C_ref;
        C_N2(j,i)=C(i,4*(p.NR+2)+j)*p.C_ref;
        C_CH4(j,i)=C(i,5*(p.NR+2)+j)*p.C_ref;
        rho(j,i)=C_H2(j,i)*p.MW(1)+C_CO(j,i)*p.MW(2)+...
            C_CO2(j,i)*p.MW(3)+C_H2O(j,i)*p.MW(4)+...
            C_N2(j,i)*p.MW(5)+C_CH4(j,i)*p.MW(6);
        T(j,i)=C(i,6*(p.NR+2)+j)*p.T_ref;
        C_tot(j,i)=C_H2(j,i)+C_CO(j,i)+...
            C_CO2(j,i)+C_H2O(j,i)+...
            C_N2(j,i)+C_CH4(j,i);
    end
end
for i=1:p.NR+1
mp(i)=((rho(i+1,end)+rho(i,end))/2)*(p.v_ref*(p.v(i+1)+p.v(i))/2)*pi()*((xr(i+1)*p.Rr)^2-(xr(i)*p.Rr)^2);
end

m_out=sum(mp); % Mass flow at the reactor outlet
m_in=p.Q*p.C_in*p.MW*p.y_in';% Mass flow at the reactor inlet


% CO COnversion
n_CO_in=p.n_in*p.y_in(2);
for i=1:length(z)
     for j=1:p.NR+1
        n_H2v(j,i)=((C_H2(j+1,i)+C_H2(j,i))/2)*...
            (p.v_ref*(p.v(j+1)+p.v(j))/2)*pi()*...
            ((xr(j+1)*p.Rr)^2-(xr(j)*p.Rr)^2);
        
        n_COv(j,i)=((C_CO(j+1,i)+C_CO(j,i))/2)*...
            (p.v_ref*(p.v(j+1)+p.v(j))/2)*pi()*...
            ((xr(j+1)*p.Rr)^2-(xr(j)*p.Rr)^2);
        
        n_CO2v(j,i)=((C_CO2(j+1,i)+C_CO2(j,i))/2)*...
            (p.v_ref*(p.v(j+1)+p.v(j))/2)*pi()*...
            ((xr(j+1)*p.Rr)^2-(xr(j)*p.Rr)^2);
        
        n_H2Ov(j,i)=((C_H2O(j+1,i)+C_H2O(j,i))/2)*...
            (p.v_ref*(p.v(j+1)+p.v(j))/2)*pi()*...
            ((xr(j+1)*p.Rr)^2-(xr(j)*p.Rr)^2);
        
        n_N2v(j,i)=((C_N2(j+1,i)+C_N2(j,i))/2)*...
            (p.v_ref*(p.v(j+1)+p.v(j))/2)*pi()*...
            ((xr(j+1)*p.Rr)^2-(xr(j)*p.Rr)^2);
        
        n_CH4v(j,i)=((C_CH4(j+1,i)+C_CH4(j,i))/2)*...
            (p.v_ref*(p.v(j+1)+p.v(j))/2)*pi()*...
            ((xr(j+1)*p.Rr)^2-(xr(j)*p.Rr)^2);
        
        y_CO(j,i)=n_COv(j,i)/(p.n_in*(1-p.y_in(4)));
        
     end
     n_H2(i)=sum(n_H2v(:,i));
     n_CO(i)=sum(n_COv(:,i));
     n_CO2(i)=sum(n_CO2v(:,i));
     n_H2O(i)=sum(n_H2Ov(:,i));
     n_N2(i)=sum(n_N2v(:,i));
     n_CH4(i)=sum(n_CH4v(:,i));
     
     n_total(i)=n_H2(i)+n_CO(i)+n_CO2(i)+n_H2O(i)+n_N2(i)+n_CH4(i);
     
     % CO conversion
     X_CO(i)=(n_CO_in-n_CO(i))/n_CO_in;
     % CO dry gas concentration
     y_COdry(i)=n_CO(i)/(n_total(i)-n_H2O(i));
end
% Outlet CO molar fraction in the dry gas
y_CO_dry=y_COdry(end);

% Outlet molar fraction
y_out=[n_H2(end) n_CO(end) n_CO2(end) n_H2O(end) n_N2(end) n_CH4(end)]/n_total(end);

%% Inert layer(s)
inert=zeros(1,length(z_in));
for i=1:length(z_in)
    if length(z_out)<length(z_in)
        z_out(length(z_in))=z(end);
    end
    inert(i)=z_out(i)*p.L-z_in(i)*p.L;
end
L_inert=sum(inert);

L=z(end)*p.L; % Tube total length (m)

P_out=p.P_in+(L*p.dPdz)/1e5; % bar 
xz=z; % Dimentionless length (-)

% Heat flux dependence of the axial coordinate (W/m^2)
qheat=C(:,end); % Total heat (W)
for i=1:(length(z)-1)
    q_flux(i)=(qheat(i+1)-qheat(i))/((z(i+1)*p.L-z(i)*p.L)*pi()*p.dt);
end

% Figures
figure(1)
subplot(2,2,1)
plot(xz,X_CO,'k')
ylabel('CO conversion')
xlabel('Dimentionless Tube Length')
set(gca,'Ytick',[0.3 0.6 0.9 1])
subplot(2,2,[2,4])
surf(xr,xz,T')
xlabel('Dimentionless Tube Radius')
ylabel('Dimentionless Tube Length')
shading interp
view(2)
h1 = colorbar;
h1.Label.String = 'Temperature (K)';
subplot(2,2,3)
plot(xz(1:end-1),q_flux,'k')
xlabel('Reactor dimentionless length')
ylabel('Heat Flux (W/m^2)')
end

function dVar=diff(z,x0,p,xr,AR,BR)
global check count z_in z_out 
C_H2=x0(1:(p.NR+2),1);
C_CO=x0((p.NR+2)+1:2*(p.NR+2),1);
C_CO2=x0(2*(p.NR+2)+1:3*(p.NR+2),1);
C_H2O=x0(3*(p.NR+2)+1:4*(p.NR+2),1);
C_N2=x0(4*(p.NR+2)+1:5*(p.NR+2),1);
C_CH4=x0(5*(p.NR+2)+1:6*(p.NR+2),1);
T=x0(6*(p.NR+2)+1:7*(p.NR+2),1);

Pz=(p.P_in+z*p.P_drop); % bar
Fpz=(Pz/1.01325)^(0.5-(Pz/1.01325)/500); % Pressure correction factor

q=p.U*(p.T_ref*T(end)-p.Tcf); % Heat flux (W/m2)

for i=1:(p.NR+2)

    p.Ct(i)=C_H2(i)+C_CO(i)+C_CO2(i)+C_H2O(i)+C_N2(i)+C_CH4(i);

    p.y(:,i)=[C_H2(i);C_CO(i);C_CO2(i);C_H2O(i);C_N2(i);C_CH4(i)]./p.Ct(i);
    
    rho=p.C_ref*p.Ct(i)*p.MW*p.y(:,i); % (kg/m3)

    p.dC(i)=(dCdT(T(i)*p.T_ref,Pz,p.y(:,i)'))*1000; % mol/m3/K

    [~,p.Dmi(:,i)]=diffusion_c(p.T_ref*T(i),Pz,p.y(:,i),1); % m2/s
    
    [~,p.mu(i)]=viscosity_mixture(T(i)*p.T_ref,Pz,p.y(:,i)'); % (Pa.s)

    p.Cp(i)=cp_real(T(i)*p.T_ref,Pz,p.y(:,i)'); %(J/mol/K)

    [~,p.lambF(i)]=thermal_conductivity_mix(T(i)*p.T_ref,Pz,p.y(:,i)'); % (W/m/K)

    p.Re(i)=rho*p.dp*p.v(i)*p.v_ref/p.mu(i);
    
    p.Pr(i)=p.Cp(i)*(1/(p.MW*p.y(:,i)))*p.mu(i)/p.lambF(i);
    
    p.lamb0(i)=p.lambF(i)*(p.lamb_cat/p.lambF(i))^(0.28-0.757*log10(p.eps(i))-...
        0.057*log10(p.lamb_cat/p.lambF(i))); % (W/m/K)
    
    p.lambDr(i)=p.lambF(i)*p.Re(i)*p.Pr(i)/(8.65*(1+19.4*((p.dp/p.dt)^2))); % (W/m/K)
    
    p.lambR(i)=p.lamb0(i)+p.lambDr(i); % (W/m/K)

    p.Dr(:,i)=((1-sqrt(1-p.eps(i))).*p.Dmi(:,i)+p.v(i)*p.v_ref*p.dp/8)/p.D_ref; % m2/s
             
    p.mDH(i)=heatofreac(T(i)*p.T_ref,Pz,p.y(:,i));
    
    p.K(i)=equiconst(T(i)*p.T_ref,Pz,p.y(:,i));
  
    p.kr(i)=5904*(1-0.55)*2.96e5*exp(-47400/(p.Rg*(p.T_ref*T(i))));
    
    if p.dp==2e-3
        p.eta(i)=p.dp1(1)+p.dp1(2)*p.y(2,i)+p.dp1(3)*p.y(4,i)...
        +p.dp1(4)*p.T_ref*T(i)+p.dp1(5)*(p.T_ref*T(i))^2+...
        p.dp1(6)*p.v_ref*p.v(i)+p.dp1(7)*p.T_ref*T(i)*p.y(2,i)...
        +p.dp1(8)*p.T_ref*T(i)*p.y(4,i);
    elseif p.dp==4e-3
        p.eta(i)=p.dp2(1)+p.dp2(2)*p.y(2,i)+p.dp2(3)*p.y(4,i)...
        +p.dp2(4)*p.T_ref*T(i)+p.dp2(5)*(p.T_ref*T(i))^2+...
        p.dp2(6)*p.v_ref*p.v(i)+p.dp2(7)*p.T_ref*T(i)*p.y(2,i)...
        +p.dp2(8)*p.T_ref*T(i)*p.y(4,i);
    elseif p.dp==6e-3
        p.eta(i)=p.dp3(1)+p.dp3(2)*p.y(2,i)+p.dp3(3)*p.y(4,i)...
        +p.dp3(4)*p.T_ref*T(i)+p.dp3(5)*(p.T_ref*T(i))^2+...
        p.dp3(6)*p.v_ref*p.v(i)+p.dp3(7)*p.T_ref*T(i)*p.y(2,i)...
        +p.dp3(8)*p.T_ref*T(i)*p.y(4,i);
    elseif p.dp==8e-3
        p.eta(i)=p.dp4(1)+p.dp4(2)*p.y(2,i)+p.dp4(3)*p.y(4,i)...
        +p.dp4(4)*p.T_ref*T(i)+p.dp4(5)*(p.T_ref*T(i))^2+...
        p.dp4(6)*p.v_ref*p.v(i)+p.dp4(7)*p.T_ref*T(i)*p.y(2,i)...
        +p.dp4(8)*p.T_ref*T(i)*p.y(4,i);
    elseif p.dp==10e-3
        p.eta(i)=p.dp5(1)+p.dp5(2)*p.y(2,i)+p.dp5(3)*p.y(4,i)...
        +p.dp5(4)*p.T_ref*T(i)+p.dp5(5)*(p.T_ref*T(i))^2+...
        p.dp5(6)*p.v_ref*p.v(i)+p.dp5(7)*p.T_ref*T(i)*p.y(2,i)...
        +p.dp5(8)*p.T_ref*T(i)*p.y(4,i);
    elseif (p.dp > 2e-3) && (p.dp < 4e-3)
        eta1=p.dp1(1)+p.dp1(2)*p.y(2,i)+p.dp1(3)*p.y(4,i)...
        +p.dp1(4)*p.T_ref*T(i)+p.dp1(5)*(p.T_ref*T(i))^2+...
        p.dp1(6)*p.v_ref*p.v(i)+p.dp1(7)*p.T_ref*T(i)*p.y(2,i)...
        +p.dp1(8)*p.T_ref*T(i)*p.y(4,i);
        dp1=2e-3;
        
        eta2=p.dp2(1)+p.dp2(2)*p.y(2,i)+p.dp2(3)*p.y(4,i)...
        +p.dp2(4)*p.T_ref*T(i)+p.dp2(5)*(p.T_ref*T(i))^2+...
        p.dp2(6)*p.v_ref*p.v(i)+p.dp2(7)*p.T_ref*T(i)*p.y(2,i)...
        +p.dp2(8)*p.T_ref*T(i)*p.y(4,i);        
        dp2=4e-3;
    
        p.eta(i)=(p.dp-dp1)*(eta2-eta1)/(dp2-dp1)+eta1;
    elseif (p.dp > 4e-3) && (p.dp < 6e-3)
        eta1=p.dp2(1)+p.dp2(2)*p.y(2,i)+p.dp2(3)*p.y(4,i)...
        +p.dp2(4)*p.T_ref*T(i)+p.dp2(5)*(p.T_ref*T(i))^2+...
        p.dp2(6)*p.v_ref*p.v(i)+p.dp2(7)*p.T_ref*T(i)*p.y(2,i)...
        +p.dp2(8)*p.T_ref*T(i)*p.y(4,i);
        dp1=4e-3;

        eta2=p.dp3(1)+p.dp3(2)*p.y(2,i)+p.dp3(3)*p.y(4,i)...
        +p.dp3(4)*p.T_ref*T(i)+p.dp3(5)*(p.T_ref*T(i))^2+...
        p.dp3(6)*p.v_ref*p.v(i)+p.dp3(7)*p.T_ref*T(i)*p.y(2,i)...
        +p.dp3(8)*p.T_ref*T(i)*p.y(4,i);
        dp2=6e-3;
    
        p.eta(i)=(p.dp-dp1)*(eta2-eta1)/(dp2-dp1)+eta1;
    elseif (p.dp > 6e-3) && (p.dp < 8e-3)
        eta1=p.dp3(1)+p.dp3(2)*p.y(2,i)+p.dp3(3)*p.y(4,i)...
        +p.dp3(4)*p.T_ref*T(i)+p.dp3(5)*(p.T_ref*T(i))^2+...
        p.dp3(6)*p.v_ref*p.v(i)+p.dp3(7)*p.T_ref*T(i)*p.y(2,i)...
        +p.dp3(8)*p.T_ref*T(i)*p.y(4,i);
        dp1=6e-3;

        eta2=p.dp4(1)+p.dp4(2)*p.y(2,i)+p.dp4(3)*p.y(4,i)...
        +p.dp4(4)*p.T_ref*T(i)+p.dp4(5)*(p.T_ref*T(i))^2+...
        p.dp4(6)*p.v_ref*p.v(i)+p.dp4(7)*p.T_ref*T(i)*p.y(2,i)...
        +p.dp4(8)*p.T_ref*T(i)*p.y(4,i);
        dp2=8e-3;
    
        p.eta(i)=(p.dp-dp1)*(eta2-eta1)/(dp2-dp1)+eta1;
    else 
        eta1=p.dp4(1)+p.dp4(2)*p.y(2,i)+p.dp4(3)*p.y(4,i)...
        +p.dp4(4)*p.T_ref*T(i)+p.dp4(5)*(p.T_ref*T(i))^2+...
        p.dp4(6)*p.v_ref*p.v(i)+p.dp4(7)*p.T_ref*T(i)*p.y(2,i)...
        +p.dp4(8)*p.T_ref*T(i)*p.y(4,i);
        dp1=8e-3;

        eta2=p.dp5(1)+p.dp5(2)*p.y(2,i)+p.dp5(3)*p.y(4,i)...
        +p.dp5(4)*p.T_ref*T(i)+p.dp5(5)*(p.T_ref*T(i))^2+...
        p.dp5(6)*p.v_ref*p.v(i)+p.dp5(7)*p.T_ref*T(i)*p.y(2,i)...
        +p.dp5(8)*p.T_ref*T(i)*p.y(4,i);
        dp2=10e-3;
    
        p.eta(i)=(p.dp-dp1)*(eta2-eta1)/(dp2-dp1)+eta1; 
    end
    
    p.r(i)=p.eta(i)*(Fpz*((p.kr(i)*(p.y(2,i)*p.y(4,i)-p.y(3,i)*p.y(1,i)/p.K(i)))))*(1000/3600); % mol/m3.s
end

%% Setting inert layer
if (T(1)*p.T_ref>523) && (check==0)
        p.r=0*p.r;
        check=1;
        count=count+1;
        z_in(count)=z; % Recording where the inert layer begins
elseif (T(1)*p.T_ref>p.T_inert_cool) && (check==1)
          p.r=0*p.r;
         z_out(count)=z; % Recording where the inert layer ends
else
        check=0;
end
    
%%        
dC_H2dz=zeros((p.NR+2),1);
dC_COdz=zeros((p.NR+2),1);
dC_CO2dz=zeros((p.NR+2),1);
dC_H2Odz=zeros((p.NR+2),1);
dC_N2dz=zeros((p.NR+2),1);
dC_CH4dz=zeros((p.NR+2),1);
dTdz=zeros((p.NR+2),1);


%% H2
% At r = 0 
dC_H2dz(1,1)=-(p.D_ref*p.Dr(1,1)/p.Rr)*AR(1,:)*C_H2(1:(p.NR+2),1);

% At 0 < r < 1 
for j=2:p.NR+1
    dC_H2dz(j,1)=(p.L/(p.v(j)*p.v_ref))*...
        ((p.D_ref*p.Dr(1,j)/(p.Rr^2))*((1/xr(j))*AR(j,:)*C_H2(1:(p.NR+2),1)+...
        BR(j,:)*C_H2(1:(p.NR+2),1))+...
        (p.D_ref/(p.Rr^2))*(1/(xr(j)-xr(j-1)))*(p.Dr(1,j)-p.Dr(1,j-1))*(AR(j,:)*C_H2(1:(p.NR+2),1))+...
        (1-p.eps(i))*(1/p.C_ref)*p.r(j,1));
end

% At r = 1 
dC_H2dz(p.NR+2,1)=-(p.D_ref*p.Dr(1,(p.NR+2))/p.Rr)*AR(p.NR+2,:)*C_H2(1:(p.NR+2),1);

%% CO
% At r = 0 
dC_COdz(1,1)=-(p.D_ref*p.Dr(2,1)/p.Rr)*AR(1,:)*C_CO(1:(p.NR+2),1);

% At 0 < r < 1 
for j=2:p.NR+1
    dC_COdz(j,1)=(p.L/(p.v(j)*p.v_ref))*...
        ((p.D_ref*p.Dr(2,j)/(p.Rr^2))*((1/xr(j))*AR(j,:)*C_CO(1:(p.NR+2),1)+...
        BR(j,:)*C_CO(1:(p.NR+2),1))+...
        (p.D_ref/(p.Rr^2))*(1/(xr(j)-xr(j-1)))*(p.Dr(2,j)-p.Dr(2,j-1))*(AR(j,:)*C_CO(1:(p.NR+2),1))-...
        (1-p.eps(i))*(1/p.C_ref)*p.r(j,1));
end

% At r = 1 
dC_COdz(p.NR+2,1)=-(p.D_ref*p.Dr(1,(p.NR+2))/p.Rr)*AR(p.NR+2,:)*C_CO(1:(p.NR+2),1);

%% CO2
% At r = 0 
dC_CO2dz(1,1)=-(p.D_ref*p.Dr(3,1)/p.Rr)*AR(1,:)*C_CO2(1:(p.NR+2),1);

% At 0 < r < 1 
for j=2:p.NR+1
    dC_CO2dz(j,1)=(p.L/(p.v(j)*p.v_ref))*...
        ((p.D_ref*p.Dr(3,j)/(p.Rr^2))*((1/xr(j))*AR(j,:)*C_CO2(1:(p.NR+2),1)+...
        BR(j,:)*C_CO2(1:(p.NR+2),1))+...
        (p.D_ref/(p.Rr^2))*(1/(xr(j)-xr(j-1)))*(p.Dr(3,j)-p.Dr(3,j-1))*(AR(j,:)*C_CO2(1:(p.NR+2),1))+...
        (1-p.eps(i))*(1/p.C_ref)*p.r(j,1));
end

% At r = 1 
dC_CO2dz(p.NR+2,1)=-(p.D_ref*p.Dr(3,(p.NR+2))/p.Rr)*AR(p.NR+2,:)*C_CO2(1:(p.NR+2),1);

%% H2O
% At r = 0 
dC_H2Odz(1,1)=-(p.D_ref*p.Dr(4,1)/p.Rr)*AR(1,:)*C_H2O(1:(p.NR+2),1);

% At 0 < r < 1 
for j=2:p.NR+1
    dC_H2Odz(j,1)=(p.L/(p.v(j)*p.v_ref))*...
        ((p.D_ref*p.Dr(4,j)/(p.Rr^2))*((1/xr(j))*AR(j,:)*C_H2O(1:(p.NR+2),1)+...
        BR(j,:)*C_H2O(1:(p.NR+2),1))+...
        (p.D_ref/(p.Rr^2))*(1/(xr(j)-xr(j-1)))*(p.Dr(4,j)-p.Dr(4,j-1))*(AR(j,:)*C_H2O(1:(p.NR+2),1))-...
        (1-p.eps(i))*(1/p.C_ref)*p.r(j,1));
end

% At r = 1 
dC_H2Odz(p.NR+2,1)=-(p.D_ref*p.Dr(4,(p.NR+2))/p.Rr)*AR(p.NR+2,:)*C_H2O(1:(p.NR+2),1);

%% N2
% At r = 0 
dC_N2dz(1,1)=-(p.D_ref*p.Dr(5,1)/p.Rr)*AR(1,:)*C_N2(1:(p.NR+2),1);

% At 0 < r < 1 
for j=2:p.NR+1
    dC_N2dz(j,1)=(p.L/(p.v(j)*p.v_ref))*...
        ((p.D_ref*p.Dr(5,j)/(p.Rr^2))*((1/xr(j))*AR(j,:)*C_N2(1:(p.NR+2),1)+...
        BR(j,:)*C_N2(1:(p.NR+2),1))+...
        (p.D_ref/(p.Rr^2))*(1/(xr(j)-xr(j-1)))*(p.Dr(5,j)-p.Dr(5,j-1))*(AR(j,:)*C_N2(1:(p.NR+2),1)));
end

% At r = 1 
dC_N2dz(p.NR+2,1)=-(p.D_ref*p.Dr(5,(p.NR+2))/p.Rr)*AR(p.NR+2,:)*C_N2(1:(p.NR+2),1);

%% CH4
% At r = 0 
dC_CH4dz(1,1)=-(p.D_ref*p.Dr(6,1)/p.Rr)*AR(1,:)*C_CH4(1:(p.NR+2),1);

% At 0 < r < 1 
for j=2:p.NR+1
    dC_CH4dz(j,1)=(p.L/(p.v(j)*p.v_ref))*...
        ((p.D_ref*p.Dr(6,j)/(p.Rr^2))*((1/xr(j))*AR(j,:)*C_CH4(1:(p.NR+2),1)+...
        BR(j,:)*C_CH4(1:(p.NR+2),1))+...
        (p.D_ref/(p.Rr^2))*(1/(xr(j)-xr(j-1)))*(p.Dr(6,j)-p.Dr(6,j-1))*(AR(j,:)*C_CH4(1:(p.NR+2),1)));
end

% At r = 1
dC_CH4dz(p.NR+2,1)=-(p.D_ref*p.Dr(6,(p.NR+2))/p.Rr)*AR(p.NR+2,:)*C_CH4(1:(p.NR+2),1);

%% T
% At r = 0
dTdz(1,1)=-(p.lambR(1,1)*p.T_ref/p.Rr)*AR(1,:)*T(1:(p.NR+2),1);

% At 0 < r < 1 
for j=2:p.NR+1
    dTdz(j,1)=(p.L/(p.Ct(j,1)*p.C_ref*p.Cp(j,1)*p.v(j)*p.v_ref))*...
        ((1/(p.Rr^2))*(1/(xr(j)-xr(j-1)))*(p.lambR(j,1)-p.lambR(j-1,1))*(AR(j,:)*T(1:(p.NR+2),1))+...
        (p.lambR(j,1)/(p.Rr^2))*((1/xr(j))*AR(j,:)*T(1:(p.NR+2),1)+BR(j,:)*T(1:(p.NR+2),1))+...
        -p.v(j)*p.v_ref*(T(j,1)/(p.Ct(j,1)*p.C_ref))*p.dC(j,1)*(p.P_drop/p.L)+...
        (1-p.eps(i))*p.mDH(j,1)*p.r(j,1)/p.T_ref);
end

% At r = 1 
dTdz((p.NR+2),1)=(p.lambR((p.NR+2),1)*p.T_ref/p.Rr)*AR(p.NR+2,:)*T(1:(p.NR+2),1)+q;

%% Total heat exchange (W)
dqdz=p.L*q*pi()*p.dt; 

%%
dVar=[dC_H2dz
          dC_COdz
          dC_CO2dz
          dC_H2Odz
          dC_N2dz
          dC_CH4dz
          dTdz
          dqdz];  
      
progress=(z)*100;
disp(['Progress(%) = ',num2str(progress)]);

end
function [v,i,d] = events(~,x0,p,xr,~,~)

C_H2=x0(1:(p.NR+2),1);
C_CO=x0((p.NR+2)+1:2*(p.NR+2),1);
C_CO2=x0(2*(p.NR+2)+1:3*(p.NR+2),1);
C_H2O=x0(3*(p.NR+2)+1:4*(p.NR+2),1);
C_N2=x0(4*(p.NR+2)+1:5*(p.NR+2),1);
C_CH4=x0(5*(p.NR+2)+1:6*(p.NR+2),1);

 for j=1:p.NR+1
    n_H2v(j)=((C_H2(j+1)+C_H2(j))/2)*...
        (p.v_ref*(p.v(j+1)+p.v(j))/2)*pi()*...
        ((xr(j+1)*p.Rr)^2-(xr(j)*p.Rr)^2);

    n_COv(j)=((C_CO(j+1)+C_CO(j))/2)*...
        (p.v_ref*(p.v(j+1)+p.v(j))/2)*pi()*...
        ((xr(j+1)*p.Rr)^2-(xr(j)*p.Rr)^2);

    n_CO2v(j)=((C_CO2(j+1)+C_CO2(j))/2)*...
        (p.v_ref*(p.v(j+1)+p.v(j))/2)*pi()*...
        ((xr(j+1)*p.Rr)^2-(xr(j)*p.Rr)^2);

    n_H2Ov(j)=((C_H2O(j+1)+C_H2O(j))/2)*...
        (p.v_ref*(p.v(j+1)+p.v(j))/2)*pi()*...
        ((xr(j+1)*p.Rr)^2-(xr(j)*p.Rr)^2);

    n_N2v(j)=((C_N2(j+1)+C_N2(j))/2)*...
        (p.v_ref*(p.v(j+1)+p.v(j))/2)*pi()*...
        ((xr(j+1)*p.Rr)^2-(xr(j)*p.Rr)^2);

    n_CH4v(j)=((C_CH4(j+1)+C_CH4(j))/2)*...
        (p.v_ref*(p.v(j+1)+p.v(j))/2)*pi()*...
        ((xr(j+1)*p.Rr)^2-(xr(j)*p.Rr)^2);
 end
 
n_H2=sum(n_H2v);
n_CO=sum(n_COv);
n_CO2=sum(n_CO2v);
n_H2O=sum(n_H2Ov);
n_N2=sum(n_N2v);
n_CH4=sum(n_CH4v);

n_total=n_H2+n_CO+n_CO2+n_H2O+n_N2+n_CH4;

% CO dry gas concentration
y_COdry=n_CO/(n_total-n_H2O);

if y_COdry>0.01
    v = 0;
else
    v=1;
end

i = 1; % Stop the integration if v=0
d = 0; % Stop either increasing or decreasing
end 