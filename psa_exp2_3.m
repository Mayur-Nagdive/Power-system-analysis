%Y-bus formulation
% Author: Mayur Nagdive
% Date: 25/09/2021
% Singular transformation method with Mutual coupling
% -----------------------------------------------------------

clc;
clear;
close all;
%------------------------------------------------------------

disp( ' FORMATION OF BUS ADMITTANCE MATRIX\n')
nb = input ( 'enter the no of buses including the slack bus ; ' );
nl = input ( 'enter the no of lines ; ' );
mb=input ('enter the no of mutal links ; ');

for j= 1:nl
    sb(j)=input ('enter the starting bus ; ');
    eb(j)=input ('enter the ending bus ; ');
    z(j) =input ('enter the input values ; ');
    disp('--------------end of this line-------------');
end
y=zeros(nl,nb);
for i=1:nl
    k=sb(i);
    m=eb(i);
    y(i,k)=1;
    y(i,m)=-1;
end
y(:, 1) = [];
u=zeros(nl,nl);
for i=1:nl
    u(i,i)=1/z(i);
end
for j=1:mb
    e=input ('enter the first line ; ');
    f=input ('enter the second line ; ');
    g=input ('enter the value ; ');
    disp('------------end of this line----------');
       u(e,e)= z(f)/(z(e)*z(f)-g^2);
       u(f,f)=z(e)/(z(f)*z(e)-g^2); 
       u(f,e)=-g/(z(f)*z(e)-g^2);
       u(e,f)=u(f,e);
end
t=y';
x=t*u;
x;
display(x);
%end-------------------------------------------------------