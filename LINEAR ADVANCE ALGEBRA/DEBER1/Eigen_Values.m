%introdución de matrices 3x3
%A=[1 2 3;4 5 6;7 8 9
%u=[1;2;3]
%d=A(;,1)la primera fila de A
%d=A(1,;)la primera columna de A
%Comprobación de dimenciones
% n el numero de iteraciones
% u el vector inicial
function [Ma,v]=Eigen_Values(A,u,n)
B=A;
e=u;
j=length(B);%length es la dimensión más grande del arreglo
%Creación de las matrices de salida
Ma=zeros(j,n+1);
v=zeros(j,n+1);
%Creación del Vo
vo=(u)/(norm(u));
%Guardar los Vo
v(:,1)=vo;
Ma(:,1)=A*vo;

if n>0
    for k=1:n
        v(:,k+1)=((A*(v(:,k)))/(norm(A*(v(:,k)))));
        Ma(:,k+1)=A*v(:,k+1);
    end
end