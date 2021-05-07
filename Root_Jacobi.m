function r=Root_Jacobi(n,alpha,betha)
% It returns the roots of the jacobi polynomial plus the end points
r=zeros(1,n);
j=1;
x=0;
erro=1;
while erro>1e-9

    [P,dP]=Jacobi(x,n,alpha,betha);

    D=P/dP;
    y=x - D;
    erro=abs(y-x);
    j=j+1;
    x=y;

end
r(1)=x;
clear j x 

for i=2:n
    erro=1;
    x(1)=r(i-1)+1e-3;
    j=1;
    while erro>1e-8

        [P,dP]=Jacobi(x,n,alpha,betha);
         D=P/dP;
        
        for l=1:i-1
            soma(l)=1/(x-r(l));
        end
        s=sum(soma);
        y=x - D*(1/(1-D*s));
        erro=abs(y-x);
        j=j+1;
        x=y;
    end
    r(i)=x;
    clear j x
    
end
r=[0 r 1];
end