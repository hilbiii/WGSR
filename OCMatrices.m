function [x,A,B] = OCMatrices(n,alpha,betha)
% The function [A,B] = OCMatrices(n,alpha,betha)
% returns the matrices A and B and the vector x
% n is the order of the Jacobi polynomial (0<x<1)
% A is the first order derivatives matrix
% B is the second order derivatives matrix
% x are the roots of the polynomial and the end points
% For Legendre: alpha=betha=0

x=Root_Jacobi(n,alpha,betha);
for i=1:n+2
    p=1;
    v(i)=0;
    for j=1:n+2
        q=x(i)-x(j);
        v(i)=p+q*v(i);
        p=q*p;
    end
end
A=zeros(n+2,n+2);
for i=1:n+2
    for j=1:n+2
         if i==j
         else
            A(i,j)=v(i)/((x(i)-x(j))*v(j));
            A(i,i)=A(i,i)-A(i,j);
        end
    end
end

B=A^2;

end

