f= @(X)
a=1
b=n-1

n= %numero de intervalos
h=(b-a)./n;

%boole's rule
xb=linspace(a,b,n/4);
for i=1:n
    yb(i)=(2.*h./45)*(7.*f(xb(i))+32.*f(xb(i./4))+12.*f(xb(i./2))+32.*f(xb(3.*i./4))+7.*f(xb(i.*4)));
end