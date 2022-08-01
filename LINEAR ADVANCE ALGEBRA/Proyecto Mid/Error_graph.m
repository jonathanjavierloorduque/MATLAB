function Error_graph(a,b,c,d,k,T,C,U0,U1,F,U)

u0 = str2func(U0);
u1 = str2func(U1);
f = str2func(F);
u = str2func(U);

P(15) = 0;
Error(15) = 0;

for i=1:15
    P(i) = i*10;
end

for l=1:15
    hx = (b-a)/(P(l)-1);
    hy = (d-c)/(P(l)-1);
    ht = T/(k-1);
    
    uij(P(l),P(l)) = 0;
    uij_m_one(P(l),P(l)) = 0;
    uij_m_two(P(l),P(l)) = 0;

    for j=1:P(l)
        for i=P(l):-1:1
            uij_m_two(i,j) = u0(a+hy*(i-1), c+hx*(j-1));
        end
    end
    
    for j=2:P(l)-1
       for i=P(l)-1:-1:2
           uij_m_one(i,j) = (ht*u1(hy*i,hx*j))+ uij_m_two(i,j);
       end
    end
    
    for z=2:k
       for j=2:P(l)-1
           for i=P(l)-1:-1:2
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
    
    uij_syn(P(l),P(l))=0;
    
    for j=1:P(l)
        for i=P(l):-1:1
            uij_syn(i,j) = u(hy*(i-1),hx*(j-1),T);
        end
    end
    Error(l) = sqrt(immse(uij,uij_syn));
end    
    figure;
    plot(P,Error);
    title('Root mean square error');
    xlabel('Partitions');
    ylabel('Error');
    
end

