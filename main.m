function F=main
% Main file
%% Reactor operational condition
T=460;          % Temprature (K)
P=20;           % Pressure (atm)
v=0.15;         % Gas velocity (m/s)
dp=0.002;       % catalyst particle diameter (m)
dt=0.05;        % Reactor interanl tubes diamenter (m)
T_inert=475;    % Inert layer minimum temperaure (K)
H2O_CO=2;       % Water-carbon dioxide molar ratio
n_dry_in=1000;  % Feed inlet molar flow rate (mol/s)
y_dry_in=[0.34 0.47 0.13 0 0 0.06]; % Feed dry inlet composition (mol.frac.)

n_H2O_in=n_dry_in*y_dry_in(2)*H2O_CO; % Water inlet molar flow rate (mol/s)
m_H2O_in=n_H2O_in*18.015/1000; % Water inlet molar flow rate (kg/s)

n_in=n_dry_in+n_H2O_in; % Total inlet molar flow rate (mol/s)

y_H2_in=n_dry_in*y_dry_in(1)/n_in; % H2 inlet molar fraction
y_CO_in=n_dry_in*y_dry_in(2)/n_in; % CO inlet molar fraction
y_CO2_in=n_dry_in*y_dry_in(3)/n_in; % CO2 inlet molar fraction
y_H2O_in=n_H2O_in/n_in; % H2O inlet molar fraction
y_N2_in=n_dry_in*y_dry_in(5)/n_in; % N2 inlet molar fraction
y_CH4_in=n_dry_in*y_dry_in(6)/n_in; % CH4 inlet molar fraction

y_in=[y_H2_in y_CO_in y_CO2_in y_H2O_in y_N2_in y_CH4_in];

P_H2O_inlet=P*y_H2O_in; % Water pressure at the reactor inlet
P_H2O_vap=exp(12.315-4386.4/(T-14.9))-1; % Water vapor pressure (bar)
H2O_Constraint=P_H2O_inlet-P_H2O_vap; % If positive, water condenses
if H2O_Constraint>0
    disp('Pressure is too high. Water may condense inside the reactor')
    disp('Decrease inlet pressure and try again')
else
    [V,~]=PRSV(T,P,y_in);  % PRSV-EoS
    C_in=(1/V)*1000;     % Total fluid concentration at inlet (mol/m3)

    At=pi()*(dt^2)/4;   % Reactor tube cross sectional area (m2)

    n_in_tube=v*C_in*At; % Inlet tube total molar flow rate (mol/s)

    [P_out,y_out,L_inert,L]=reactor(T,P,n_in_tube,y_in,dp,dt,T_inert);

    Nt=n_in/n_in_tube; % Number of tubes

    A=pi()*dt*L*Nt;     % Total heat exchange area (m2)

    P_shell=exp(12.315-4386.4/(T-14.9))-1; % Reactor shell pressure (bar)

    %% Compression work prior to reactor inlet 
    P1=2; % Compressor inlet pressure (bar)
    % Compression ratio
    CR=P/P1;
    if CR <= 3
        if P<2.1
                Wfi=0; % Fluid power (kW)
                Wsi=0; % Shaft power (kW)
                We=0; % Eletric power (kW)
                NS=1; % Number of compression stages
                P2=[P1 P]; % Pressure Inlet/Outlet (bar)
        else
        NS=1; % Number of compression stages
        P2=[P1 P]; % Pressure Inlet/Outlet (bar)
        [~,DHf]=compression_hentalpy(P2(1),P2(2),y_dry_in);
        Wfi=(DHf*n_dry_in)/1000; % Fluid power (kW)
        Wsi=Wfi/0.78; % Shaft power (kW)
        We=Wsi/0.98; % Eletric power (kW)
        end
    elseif CR <= 9
        NS=2;
        CRi=CR^(1/NS); % Compression ratio increase by each stage
        P2=[P1 (P1*CRi) (P1*CRi*CRi)]; % Pressure Inlet/Outlet (bar)
        for i=1:(length(P2)-1)
            [T2,DHf]=compression_hentalpy(P2(i),P2(i+1),y_dry_in);
            Wfi(i)=(DHf*n_dry_in)/1000; % Fluid power (kW)
            Wsi(i)=Wfi(i)/0.78; % Shaft power (kW)
            Wei(i)=Wsi(i)/0.98; % Eletric power (kW)
            Wf=sum(Wfi); % Total fluid power (kW)
            Ws=sum(Wsi); % Total shaft power (kW)
            We=sum(Wei); % Total Eletric power (kW) 

            if i>1
                [A_HE(i-1),Q_H2O(i-1)]=heat_exchanger(T2,P2(i),y_dry_in,n_dry_in);
            end
        end
        Q_water=sum(Q_H2O);
    elseif CR <= 27
        NS=3;
        CRi=CR^(1/NS); % Increase in Pressure by each stage
        P2=[P1 (P1*CRi) (P1*CRi*CRi) (P1*CRi*CRi*CRi)]; % Pressure 
        % Inlet/Outlet (bar)
        for i=1:(length(P2)-1)
            [T2,DHf]=compression_hentalpy(P2(i),P2(i+1),y_dry_in);
            Wfi(i)=(DHf*n_dry_in)/1000; % Fluid power (kW)
            Wsi(i)=Wfi(i)/0.78; % Shaft power (kW)
            Wei(i)=Wsi(i)/0.98; % Eletric power (kW)
            Wf=sum(Wfi); % Total fluid power (kW)
            Ws=sum(Wsi); % Total shaft power (kW)
            We=sum(Wei); % Total Eletric power (kW) 

            if i>1
                [A_HE(i-1),Q_H2O(i-1)]=heat_exchanger(T2,P2(i),y_dry_in,n_dry_in);
            end
        end
        Q_water=sum(Q_H2O);
    else
        NS=4;
        CRi=CR^(1/NS); % Increase in Pressure by each stage
        P2=[P1 (P1*CRi) (P1*CRi*CRi) (P1*CRi*CRi*CRi) (P1*CRi*CRi*CRi*CRi)];
        % Pressure Inlet/Outlet (bar)
        for i=1:(length(P2)-1)
            [T2,DHf]=compression_hentalpy(P2(i),P2(i+1),y_dry_in);
            Wfi(i)=(DHf*n_dry_in)/1000; % Fluid power (kW)
            Wsi(i)=Wfi(i)/0.78; % Shaft power (kW)
            Wei(i)=Wsi(i)/0.98; % Eletric power (kW)
            Wf=sum(Wfi); % Total fluid power (kW)
            Ws=sum(Wsi); % Total shaft power (kW)
            We=sum(Wei); % Total Eletric power (kW)

            if i>1
                [A_HE(i-1),Q_H2O(i-1)]=heat_exchanger(T2,P2(i),y_dry_in,n_dry_in);
            end
        end
        Q_water=sum(Q_H2O);
    end

    %% Compression work penalty due to pressure drop inside the reactor
    y_H2_dry_out=y_out(1)/(1-y_out(4));
    y_CO_dry_out=y_out(2)/(1-y_out(4));
    y_CO2_dry_out=y_out(3)/(1-y_out(4));
    y_H2O_dry_out=0;
    y_N2_dry_out=y_out(5)/(1-y_out(4));
    y_CH4_dry_out=y_out(6)/(1-y_out(4));
    y_dry_out=[y_H2_dry_out y_CO_dry_out y_CO2_dry_out...
        y_H2O_dry_out y_N2_dry_out y_CH4_dry_out];

    if P_out < 0
        P_out=0.01;
        DHp=compression_hentalpy(P_out,P,y_dry_out);
    else
        DHp=compression_hentalpy(P_out,P,y_dry_out);
    end
    Wfp=(DHp*n_dry_in)/1000; % Fluid power (kW)
    Wsp=Wfp/0.78; % Shaft power (kW)
    Wep=Wsp/0.98; % Electric power(kW)

    %% Catalyst mass and inert mass
    m_catalyst=1441.62*At*Nt*(L-L_inert); % kg
    m_inert=1361.57*At*Nt*L_inert; % kg

    %% Economic analysis
    % Bare module cost - HE between compressors
    if NS>1
        for i=1:(NS-1)
            if P2(i+1)<5
                Fp=1;
            else
                Fp=10^(0.03881-0.11272*log10(P2(i+1))+...
                    0.08183*((log10(P2(i+1)))^2));
            end
            C_HEX(i)=(1.63+2.988*Fp)*10^(4.8306-0.85009*log10(A_HE(i))+...
                0.3187*((log10(A_HE(i)))^2));
        end
    else
        C_HEX=0;
    end
    C_HE=sum(C_HEX);

    % Bare module cost - Reactor
    P_design=max([(abs(P_shell-P2(end))) P_shell]);
    if P_design<5
        Fp=1;
    else
        Fp=10^(0.03881-0.11272*log10(P_design)+0.08183*((log10(P_design))^2));
    end
    C_reac=(1.63+2.988*Fp)*10^(4.8306-0.85009*log10(A)+0.3187*((log10(A))^2));

    for i=1:NS
        if Wfi>0
            % Bare module cost - Compressor
            C_comp(i)=2.8*10^(2.2897+1.3604*log10(Wfi(i))-...
                0.1027*((log10(Wfi(i)))^2));
            % Bare module cost - Compressor Driver
            C_driv(i)=1.5*10^(1.956+1.7142*log10(Wsi(i))-...
                0.2286*((log10(Wsi(i)))^2));
        else
            % Bare module cost - Compressor
            C_comp(i)=0;
            % Bare module cost - Compressor Driver
            C_driv(i)=0;
        end
    end
    C_compr=sum(C_comp);
    C_driver=sum(C_driv);

    % Total bare module cost
    C_BM=1.51461*(C_reac+C_compr+C_driver+C_HE);

    % Operating Costs
    % Steam cost
    if m_H2O_in < 40
        C_STEAM=(2.3e-5*m_H2O_in^(-0.9)*601.3+0.0034*P_shell^(0.05)*3.19)*...
        m_H2O_in*24*3600*330*20; 
    else
        C_STEAM=(2.3e-5*40^(-0.9)*601.3+0.0034*P_shell^(0.05)*3.19)*...
        m_H2O_in*24*3600*330*20; 
    end

    % Cooling water
    if NS>1
        if Q_water < 10
            C_CW=((7e-5+2.5e-5*Q_water^(-1))*601.3+0.003*3.19)*...
            Q_water*24*3600*330*20;
        else
            C_CW=((7e-5+2.5e-5*10^(-1))*601.3+0.003*3.19)*...
            Q_water*24*3600*330*20;
        end
    else
        C_CW=0;
    end

    % Total Operating cost
    C_O=10628.64*(We+Wep)+123.55*m_catalyst+15.48*m_inert+C_STEAM+C_CW;

    %% Total cost (US$)
    F=C_O+C_BM;
end
format long
disp(F)
format short
end