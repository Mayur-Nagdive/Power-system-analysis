% Name: Mayur Nagdive
%Sub date: 26-10-21
close all; clear; clc;

%---------Fast Decoupled Newton Raphson------%
st = clock;
b=busdat;   l=linedat;
nb=length(b(:,1));
n1=find(b(:,2)>1);  % finds index of non slack buses
n2=find(b(:,2)>2);  % finds index of pq buses
Y=ybusf(l);         %
B1 =imag(Y(n1,n1));
B11=imag(Y(n2,n2));
Ps=b(:,5)-b(:,7);   dP=Ps;  % schedule powers and
Qs=b(:,6)-b(:,8);   dQ=Qs;  % initial values
d=b(:,4);V=b(:,3);
while (max(abs(dP(n1)))>0.0001||max(abs(dQ(n2)))>0.0001)
    S=(V.*exp(-1i*d)).*(Y*(V.*exp(1i*d)));  % S=(V-d)x(Yx(V?d))
    dP=Ps-real(S);  dQ=Qs+imag(S);          % S=V*I || S=P-jQ
    d(n1)=d(n1)-B1\ dP(n1) ./V(n1);
    V(n2)=V(n2)-B11\dQ(n2)./V(n2);
end
ste = clock;

% -----------------------Diplaying Outputs---------------------%
disp('                     Fast Decoupled Load Flow Analysis');
disp('                     Report of power flow calculations')
disp(date);
disp('The  Y-bus is:')
disp(Y)
fprintf('Solution time              : %g sec.\n',etime(ste,st))
display(V);
D=d*180/pi;
display(D);

%--------------------Functions-------------------------%
function d = busdat
Slack=1;    PV=2;   PQ=3;
%     |Bus | Type | Vsp |del| PGi | QGi | PLi | QLi |
d = [   1   Slack  1.06   0    0     0     0       0;
        2    PV    1.0    0   .40   .30   .20    .10;
        3    PQ    1.0    0    0     0    .45    .15;
        4    PQ    1.0    0    0     0    .40    .05;];
end

function d = linedat
%     |FromBus|ToBus|Impedance|LineCharging|
d = [    1      2    .02+.06j     .030j;
         1      3    .08+.24j     .025j;
         2      3    .06+.18j     .020j;
         2      4    .06+.18j     .020j;
         3      4    .01+.03j     .010j;];
end

function Y = ybusf(l)
fb = l(:,1);    tb = l(:,2);
y = 1./l(:,3);   b = l(:,4) + y;
nb = max(max(fb,tb));	Y = zeros(nb);
Y(nb*(fb-1)+tb)=-y; Y(nb*(tb-1)+fb)=-y;
    for m=1:nb
        Y(m,m)=Y(m,m)+sum(b(fb==m))+sum(b(tb==m));
    end
end
%-----------------------end------------------------------%