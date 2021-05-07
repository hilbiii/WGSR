function [v,dPdz]=momentum(T,P,y,Q,dp,dr,N)
% Function 'momentum' calculates the velocity profile (v) in m/s as a
% function of the radial position of the recator and the pressure drop in
% Pa/m (dP/dz). Since the momentum balance model considers that phisical
% properties variation cause no effect on the velocity field, velocity is 
% function of the radial position due to porosity variation only and the
% pressure is a function of axial position only.
% Therefore, the velocity profile is found by guessing the
% pressure drop until the volumetric flow rate inside the packed bed is
% equal to the inlet volumetric flow rate (actually, the difference is less
% than 1e-8). This function requires 5 inputs: Temperature in K, pressure
% in bar, molar fraction, inlet volumetric flow rate m3/s, catalyst 
% particle diamter in m, reactor tube diameter in m and number of internal
% collocation points in the radial coordinate.
% [v,dPdz]=momentum(400.25,1.01325,[0.435 0.08 0.11 0.319 0.006 0.05],...
% 7.55357e-4,5e-4,0.072)

% Transfering input variables
p.P_in=P;           % bar
p.T_in=T;           % K
p.y_in=y;         	% -
p.dp=dp;            % m
p.Rr=dr/2;          % m
p.Q_in=Q;           % m3
p.MW=[2.016,28.01,44.01,18.015,28.014,16.043]; % Molecular weight (kg/kmol)
p.g= 9.81;          % m/s2
p.v_ref=0.001;      % m/s

% Fluid inlet density
[Vm,~]=PRSV(p.T_in,p.P_in,p.y_in);
p.rho=p.MW*p.y_in'/Vm; % kg/m3

% Fluid inlet viscosity
[~,p.mu]=viscosity_mixture(p.T_in,p.P_in,p.y_in); % Pa.s

% Setting orthogonal collocation matrices
p.N=N; % Number of internal collocation points
[x,A,B] = OCMatrices(p.N,0,0);

% Porosity as a function of the radial coodinate
p.epsilon=zeros(1,p.N+2);
for i=1:p.N+2
    p.epsilon(i)=0.4*(1+1.36*exp(-5*(p.Rr-x(i)*p.Rr)/p.dp));
end

% Setting velocity profile intial guesses
X0=0.001*ones(1,p.N+2);
% load X0
% Setting solver options
op=optimoptions('lsqnonlin','TolFun',1e-10,...
'TolX',1e-10,'MaxFunEvals',6e4,...
'MaxIter',6e4,'Display','off');

% Guessing the pressure drop dPdz=K
K(1)=-8e+03;% 18 pt
Q=0;
residual=abs(Q-p.Q_in);
i=0;
Inc=1;
while residual>1e-8
    i=i+1;
    v = lsqnonlin(@bal_momentum,X0,zeros(1,p.N+2),[],op,p,x,A,B,K(i));
    % Calculating the volumetric flow rate for the velocity profile 
    % calculated with the guessed dPdz 
    Q=zeros(1,p.N+1);
    for j=1:p.N+1
        Q(j)=(p.v_ref*(v(j+1)+v(j))/2)*pi()*...
            ((x(j+1)*p.Rr)^2-(x(j)*p.Rr)^2);
    end
    Q=sum(Q);
    residual=abs(Q-p.Q_in);
    % Checking the difference between the calculated Q and inlet Q
    if Q>p.Q_in
        K(i+1)=K(i)+100/Inc;% Updating dPdz
    else
        K(i+1)=K(i)-100/Inc;% Updating dPdz
            if i>3
                if (K(i)-K(i-2)==0)
                    Inc=Inc*10;
                    K(i+1)=K(i)-100/Inc;% Updating dPdz
                end
            end
    end
    X0=v; % Updating solver intial guess
end
v=p.v_ref*v;
dPdz=K(end);
end

function F=bal_momentum(X0,p,x,A,B,K)
v=X0;
F=zeros(1,p.N+2);
% At x=0
F(1)=A(1,:)*v';
% At 0 < x < 1
for i=2:p.N+1
    F(i)=-K+p.rho*p.g+(p.mu*p.v_ref/(p.epsilon(i)*...
        x(i)*(p.Rr^2)))*A(i,:)*v'...
        +(p.mu*p.v_ref/(p.epsilon(i)*(p.Rr^2)))*B(i,:)*v'...
        -150*((1-p.epsilon(i))^2)*p.mu*p.v_ref*v(i)/...
        ((p.epsilon(i)^3)*(p.dp^2))...
        -1.75*(1-p.epsilon(i))*((p.v_ref*v(i))^2)*...
        p.rho/(p.dp*(p.epsilon(i)^3));
end
% At x=1
F(p.N+2)=v(p.N+2); 
end


