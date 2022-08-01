function Midterm_temp(a,b,c,d,n,m,k,T,C) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f = @(x,y,t) 2*(a-x)*(x-b)*(c-y)*(y-d)-2*(t^2)*(y-c)*(y-d)-2*(t^2)*(x-a)*(x-b);
u0 = @(x,y) 0;
u1 = @(x,y) 0;
u = @(x,y,t) (x-a)*(y-c)*(x-b)*(y-d)*(t^2);

uij(m,n) = 0;
uij_m_one(m,n) = 0;
uij_m_two(m,n) = 0;
uij_syn(m,n) = 0;
Error(m,n) = 0;

Pj = linspace(a,b,n);
Pi = linspace(c,d,m);

hx = (b-a)/(n-1);
hy = (d-c)/(m-1);
ht = T/(k-1);

    for j=1:n
        for i=m:-1:1
            uij_m_two(i,j) = u0(a+hy*(i-1), c+hx*(j-1));
        end
    end
    
    for j=2:n-1
       for i=m-1:-1:2
           uij_m_one(i,j) = (ht*u1(hy*i,hx*j))+ uij_m_two(i,j);
       end
    end
    
    for z=2:k
       for j=2:n-1
           for i=m-1:-1:2
               dux = (uij_m_one(i,j-1)-2*uij_m_one(i,j)...
               +uij_m_one(i,j+1))/(hx^2);
           
               duy = (uij_m_one(i-1,j)-2*uij_m_one(i,j)...
               +uij_m_one(i+1,j))/(hy^2);
               
               dut = uij_m_two(i,j)-2*uij_m_one(i,j);
               
               uij(i,j) = (C^2)*(ht^2)*(f(hy*i,hx*j,ht*(z-1))+dux+duy)-dut;
           end
       end
       uij_m_two = uij_m_one;
       uij_m_one = uij;
    end
    
    for j=1:n
        for i=m:-1:1
            uij_syn(i,j) = u(hy*(i-1),hx*(j-1),T);
        end
    end
    
    for j=1:n
        for i=m:-1:1
            Error(i,j) = abs(uij_syn(i,j)-uij(i,j));
        end
    end
    figure;
    plot(Error);
    disp(Error);
    
    figure;
    surf(Pj,Pi,uij);
    figure;
    surf(Pj,Pi,uij_syn);
end
