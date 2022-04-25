
%Modelling: parameters of the transmission line
%Author: Mayur Nagdive
%Date: 20/08/2021
%---------------------------------------------------------------------------

close all;
clear;
clc;
%----------------------------------------------------------------------------

n=input('Enter the number of the buses present in the system \n');
l=input("enter the length of the conductor \n");
d=input("enter the diameter of the conductor \n");
f=input('enter the frequncy of the system \n');
rho=input('enter the resistivity of the  conductor \n');
a= input('Spacing between wire(1): \n');
b= input('Spacing between wire(2): \n');
c= input('Spacing between wire(3): \n');
%------------------------------------------------------------------------------

r=d/2;
g=(a*b*c)^(1/3);
q=log(g/r)/(log(exp(1)));
R= ((rho)/(pi*(r^2)))*l;
L= 2*((q)+0.25)*l*(10^(-7));
e=8.854187817*(10^(-12));
C=((2*pi*e)/(q))*l;
X=(2*pi*50*L);
Y=1i*(2*pi*50*C);
Z=(R+1i*X);

%--------------------------------------------------------------------------------

fprintf('Resistance (in ohm) = %f \n' ,R );
fprintf('Inducance (in henry)= %f \n' , L);
fprintf('Capacitance (in Farad) = %f \n' , C);
% ------------------------------------------------------------------------------

if l<8000
    fprintf('Its a short Transmission line:\n');
    fprintf('The ABCD PARAMENTERS ARE:\n');
    A=1;
    B=(R+1i*X);
    C=0;
    D=1;
    fprintf('A is: %f \n',A);
    fprintf('B is: %f \n',B);
    fprintf('C is: %f \n',C);
    fprintf('D is: %f \n',D);
     

elseif l<160000
    fprintf('Its a medium Transmission line:\n');
    fprintf('For nominal pi model The ABCD parametres are as follows: \n');
    A=1+(Z*Y)/2;
    B=Z;
    C=Y*(1+(Y*Z/4));
    D=1+(Z*Y)/2;
    fprintf('A is: %f \n',A);
    fprintf('B is: %f \n',B);
    fprintf('C is: %f \n',C);
    fprintf('D is: %f \n',D);
     
     fprintf('For nominal T model The ABCD parametres are as follows:\n');
     A=1+(Z*Y)/2;
    B=Z*(1+(Y*Z/4));
    C=Y;
    D=1+(Z*Y)/2;
    fprintf('A is: %f \n',A);
    fprintf('B is: %f \n',B);
    fprintf('C is: %f \n',C);
    fprintf('D is: %f \n',D);
    
else
    fprintf('Its a Long Transmission line:\n')
    fprintf('The ABCD PARAMENTERS ARE:\n')
    v=sqrt(Z/Y);
    o=(sqrt(Z*Y));  %propogation constant*length
    A= cosh(o);
    B=v*(sinh(o));
    C= (sinh(o))/v;
    D=cosh(o);
    fprintf('A is: %f \n',A);
    fprintf('B is: %f \n',B);
    fprintf('C is: %f \n',C);
    fprintf('D is: %f \n',D);
    fprintf('V is: %f \n',v);
    fprintf('o is: %f',o);
end