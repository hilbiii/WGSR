function [T2,DHc]=compression_hentalpy(P1,P2,y)
% It caculates the enthalpy change of a gas mixture during compression from
% P1 to P2, assuming that the gas intial temperature is 305 K.
% y=[H2,CO,CO2,H2O,N2,CH4] is the componets molar fraction array

T1=305;
p.R=8.314; % Ideal gas constant
%% Cp mixture ideal gas
p.a0=sum([2.883 3.912 3.259 4.395 3.539 4.568].*y);
p.a1=sum([3.681 -3.913 1.356 -4.186 -0.261 -8.975]*1e-3.*y);
p.a2=sum([-0.772 1.182 1.502 1.405 0.007 3.631]*1e-5.*y);
p.a3=sum([0.692 -1.302 -2.374 -1.564 0.157 -3.407]*1e-8.*y);
p.a4=sum([-0.213 0.515 1.056 0.632 -0.099 1.091]*1e-11.*y);
x0=T1+50; % Guess outlet temperature

% Setting solver options
op=optimoptions('fsolve','TolFun',1e-15,...
'TolX',1e-15,'MaxFunEvals',10e6,...
'MaxIter',1e6,'Display','off');

T2 = fsolve(@temperature,x0,op,p,P1,P2,T1,y);

DHc=hentalpy_delta(T1,T2,P1,P2,y);
end

function res=temperature(x0,p,P1,P2,T1,y)
T2=x0;
Sig=p.R*(p.a0*log(T2/T1)+p.a1*(T2-T1)+(p.a2/2)*(T2^2-T1^2)+...
    +(p.a3/3)*(T2^3-T1^3)+(p.a4/4)*(T2^4-T1^4))-p.R*log(P2/P1); % Ideal gas
Sres1=residual_entropy(T1,P1,y); % Residual at T1 and P1
Sres2=residual_entropy(T2,P2,y); % Residual at T2 and P2

res=Sig+Sres2-Sres1;

end
