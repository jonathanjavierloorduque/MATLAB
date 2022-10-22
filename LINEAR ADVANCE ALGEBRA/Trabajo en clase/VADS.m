function [u] = VADS(c,f,l,n)


h=1/n;
A(n-1,n-1)=0;
b(n-1)=0;

w=@(x,y)1-1/h*abs(x-y*h);

b(1)=boole(@(x)f(x).*w(x,1), 0,2*h);
A(1,1)=2/h*boole(@(x)c(x)*w(x,1)^2, 0, 2*h);
A(1,2)=-1/h+boole(@)(x)c(x)*w(x,1)*w(w,2),h,2*h);

for k=2;n-2;
    b(k)=booble(@(x)f(x)*w(x,k))

end

u=A\b;

end

