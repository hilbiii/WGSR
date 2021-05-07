function [P,dP]=Jacobi(x,n,alpha,betha)
%% Algorithm used to calculate the roots of the Jacobi polynomial (x).
% P is the polynomial and dP is its derivate.
a=alpha;
b=betha;

if n==0
    P=1;
    dP=0;
else
    for i=1:n
        if i==1
            g(i)=(a+1)/(a+b+2);
            h(i)=0;
            P(i)=(x-g(i))*1;
            dP(i)=1;
        elseif i==2
            g(i)=0.5*(1-((a^2-b^2)/(((2*i+a+b-1)^2)-1)));
            h(i)=(a+1)*(b+1)/(((a+b+2)^2)*(a+b+3));
            P(i)=(x-g(i))*P(i-1)-h(i)*1;
            dP(i)=(x-g(i))*dP(i-1)-h(i)*0+P(i-1);
        else
            g(i)=0.5*(1-((a^2-b^2)/(((2*i+a+b-1)^2)-1)));
            h(i)=(i-1)*(i+a-1)*(i+b-1)*(i+a+b-1)/((2*i+a+b-1)*((2*i+a+b-2)^2)*(2*i+a+b-3));
            P(i)=(x-g(i))*P(i-1)-h(i)*P(i-2);
            dP(i)=(x-g(i))*dP(i-1)-h(i)*dP(i-2)+P(i-1);
        end
    end
    P=P(end);
    dP=dP(end);
end

end