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
%[3.741657386773941,3.738207263968532,3.738664012491395;8.552359741197580,8.460153281612993,8.466974381230510;13.363062095621220,13.182099299257452,13.195284749969629]