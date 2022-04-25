%Y-bus formulation
% Author: Mayur Nagdive
% Date: 25/09/2021
% step-by -step method
% -----------------------------------------------------------

close all;
clear;
clc
%------------------------------------------------------------

disp('----------Ybus formation code--------');
nbranch=input('enter the number of branches in system = ');
disp('enter line data');
for n=1:1:nbranch
    fb=input('enter from bus = ');
    tb=input('enter to bus = ');
    r=input('enter value of resistance = ');
    x=input('enter value of reactance = ');   
    B=input('enter the value of line charging admittance(b/2) = ');
    disp('------------end of this branch------------')
    
    %calculating admittances
    z=r+1i*x;
    y=1./z;
    Ldata(n,:)=[fb tb r x B y];
end


fb=Ldata(:,1);
tb=Ldata(:,2);
r=Ldata(:,3);
x=Ldata(:,4);
b=Ldata(:,5);
y=Ldata(:,6);

b=1i*b;

nbus=max(max(fb),max(tb));
Y=zeros(nbus,nbus);

for k=1:nbus
    bs=input('enter the bus at which ground reactor is linked = ');
    x0=input('enter value of ground reactor = ');
    disp('---------------------------')
    if x0~=0
        y0=1./(1i*x0);
    else
        y0=0;
    end      
   Yreact(k,:)=[bs y0];
end
bs=Yreact(:,1);
y0=Yreact(:,2);

% off diagonal element
for k=1:nbranch
    Y(fb(k),tb(k))=Y(fb(k),tb(k))-y(k);
    Y(tb(k),fb(k))=Y(fb(k),tb(k));
end

% diagonal element
for m=1:nbus
    for n=1:nbranch
        if fb(n)==m 
            Y(m,m)=Y(m,m)+y(n)+b(n);
        elseif tb(n)==m          
            Y(m,m)=Y(m,m)+y(n)+b(n);
        end
    end
end

for k=1:nbus
    Y(k,k)=Y(k,k)+y0(k);
end

Yb=Y;
display(Yb);
%end-------------------------------------------------------