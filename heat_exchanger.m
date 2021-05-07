function [A_HE,Q_H2O]=heat_exchanger(T_gas,P_gas,y_dry_in,n_dry_in)
% It calculates the heat exhanger(m2) area and heat duty(kW)
U=70; % W/m2.K
cp_H2O=4180; % J/kg.K
rho_H2O=1000; % kg/m3
Tc_in=298; % K
Tc_out=313; % K
Th_in=T_gas; % K
Th_out=305; % K
DT1=Th_in-Tc_out; % K
DT2=Th_out-Tc_in; % K

DTLM=(DT1-DT2)/log(DT1/DT2);
DH_gas=hentalpy_delta(Th_out,Th_in,P_gas,P_gas,y_dry_in);

Q=DH_gas*n_dry_in; % W
Q_H2O=Q/(cp_H2O*(Tc_out-Tc_in)*rho_H2O); % m3/s
A_HE=Q/(U*DTLM); % m2

end