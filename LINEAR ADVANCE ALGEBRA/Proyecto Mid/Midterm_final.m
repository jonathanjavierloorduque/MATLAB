function Midterm_final(a,b,c,d,n,m,k,T,C) 
prompt = 'Enter u0, eg: @(x,y) x+y+1 ';
u0 = str2func(input(prompt, "s"));
prompt1 = 'Enter u1, eg: @(x,y) x+y+1 ';
u1 = str2func(input(prompt1, "s"));
prompt2 = 'Enter f, eg: @(x,y,t) x+y+t+1 ';
f = str2func(input(prompt2, "s"));
prompt3 = 'Enter u, eg: @(x,y,t) x+y+t+1\nIf you do not have the u formula, press "Enter" ';

try
    u = str2func(input(prompt3, "s"));
catch
    u = '*';
end


uij(m,n) = 0;
uij_m_one(m,n) = 0;
uij_m_two(m,n) = 0;
uij_syn(m,n) = 0;
Error(m,n) = 0;
A_matrix((m-2)*(n-2),(n-2)*(m-2)) = 0;
multiples(n-2) = 0;

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
    figure;
    surf(Pj,Pi,uij);
    title('Approximated Solution');
    xlabel('X');
    ylabel('Y');
    zlabel('T');
    
    
    if isa(u,'function_handle')
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
        surf(Pj,Pi,uij_syn);
        title('Exact Solution');
        xlabel('X');
        ylabel('Y');
        zlabel('T');
         
        Error_graph(a,b,c,d,k,T,C,char(u0),char(u1),char(f),char(u));
        
    end
    
    for i=1:(m-2)
        multiples(i) = i*(n-2);
    end
    
    for i=1:(m-2)*(n-2)
        for j=1:(n-2)*(m-2)
            if i==j
                A_matrix(i,j) = 1;
                if j > 1
                    A_matrix(i,j-1) = 1;
                end
                if j < (n-2)*(m-2)
                    A_matrix(i,j+1) = 1;
                end
                if j-(n-2) >= 1
                    A_matrix(i,j-(n-2)) = 1;
                end 
                if j+n-2 <= (m-2)*(n-2)
                    A_matrix(i,j+n-2) = 1;
                end
                if ismember(j,multiples) && (j < (n-2)*(m-2)) 
                    A_matrix(i,j+1) = 0;
                end
                if ismember(j-1,multiples)
                    A_matrix(i,j-1) = 0;
                end
            end
        end
    end
    
    
    
    figure;
    spy(A_matrix);
    title('Associated Matrix');
end
