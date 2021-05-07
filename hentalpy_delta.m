function DH=hentalpy_delta(T1,T2,P1,P2,y)
% It returns the enthalpy difference for a gas mixture between 2 states

p.R=8.314; % Ideal gas constant
% Cp mixture ideal gas constants
p.a0=sum([2.883 3.912 3.259 4.395 3.539 4.568].*y);
p.a1=sum([3.681 -3.913 1.356 -4.186 -0.261 -8.975]*1e-3.*y);
p.a2=sum([-0.772 1.182 1.502 1.405 0.007 3.631]*1e-5.*y);
p.a3=sum([0.692 -1.302 -2.374 -1.564 0.157 -3.407]*1e-8.*y);
p.a4=sum([-0.213 0.515 1.056 0.632 -0.099 1.091]*1e-11.*y);

% Ideal gas enthalpy
Hig=p.R*(p.a0*(T2-T1)+(p.a1/2)*(T2^2-T1^2)+(p.a2/3)*(T2^3-T1^3)+...
    (p.a3/4)*(T2^4-T1^4)+(p.a4/4)*(T2^5-T1^5)); % J/mol

% Residual enthalpy (J/mol)
H1=residual_enthalpy(T1,P1,y);
H2=residual_enthalpy(T2,P2,y);

DH=Hig+H2-H1; %(J/mol)
end