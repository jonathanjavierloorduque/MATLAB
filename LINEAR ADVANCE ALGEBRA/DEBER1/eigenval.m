%A=[1 2 3;4 5 6;7 8 9]
%u=[1;2;3]
%n=2
function [lam] = eigenval(A,K,u)
B=A;
%e=u;
[n,m] = size(A);
for i = 1:n
    for k = 1:K
        v = A*u;
        [vmax,I] = max(abs(v));
        u = v/v(I);
        a(k) = ((u.'*A)*u)/(u.'*u);
    end
    lam(i) = a(K);
    x1 = u;
    e1 = zeros(n - i + 1,1);
    e1(1) = 1;
    %w = u + norm(u)*e1;
    %H = (eye((n -i + 1) - (2*(w*w.'))/((w.'*w))));
    %B = H*(A)*H;
    %A = B(2:end,2:end);
end
%[16.117778772951187,16.116848461521510,16.116843991376460]