% Name: Mayur Nagdive
%Sub date: 26-10-21
close all; clear; clc;


 % Information about the bus matrix
  % nd   V    Ang.    Pg      Qg      PL      QL      Gs       jBs    Type
  % (1) (2)   (3)     (4)     (5)     (6)     (7)     (8)      (9)    (10)
  % Colum 11: if the bus has shunt element =1, if it hasnt shunt element =0
  
  bus=[ 1  1.06  0.000   0.00    0.00    0.00    0.00    0.000    0.000  1  0.0;
            2  1.00  0.000   0.40    0.00    0.20    0.10    0.000    0.000  2  0.0;
            3  1.00  0.000   0.00    0.00    0.45    0.15    0.000    0.000  3  0.0;
            4  1.00  0.000   0.00    0.00    0.40    0.05    0.000    0.000  3  0.0];
 
 %Information about the line matrix
%COL 1.-  From bus
%COL 2.-  to bus
%COL 3.-  R P.U.
%COL 4.-  Xl P.U.
%COL 5.-  Xc (parallel) P.U.
%COL 6.-  Type of line: 0==Line ; 1==Transformer
%COL 7.- phase shifter angle
line=[  1  2   0.02   0.06   0.0150   0.00   0.00;
           1  3   0.08   0.24   0.0125   0.00   0.00;
           2  3   0.06   0.18   0.0100   0.00   0.00;
           2  4   0.06   0.18   0.0100   0.00   0.00;
           3  4   0.01   0.03   0.0050   0.00   0.00]; 

% Data of the power system is stored in the bus and line matrix of the
% following file
nbuses=length(bus(:,1)); % number of buses of the electric power system
V=bus(:,2); Vprev=V; % Initial bus voltages
The=bus(:,3); % Initial bus angles
% Net power (Generation - Load)
P=bus(:,4)-bus(:,6);     
Q=bus(:,5)-bus(:,7); 
% ++++++++++++++ First, compute the addmitance matrix Ybus ++++++++++++++++
[Y] = Y_admi(line,bus,nbuses); % function to get the admittance matrix 

% +++++++++++++++++++++++ Start iterative process +++++++++++++++++++++++++
tolerance=1; 
iteration=0;
st=clock; % start the iteration time clock
while (tolerance > 1e-8)
    for k=2:nbuses
        PYV=0;
        for i=1:nbuses
            if k ~= i
                PYV = PYV + Y(k,i)* V(i);  % Vk * Yik
            end
        end
        if bus(k,10)==2 % PV bus
            % Estimate Qi at each iteration for the PV buses
            Q(k)=-imag(conj(V(k))*(PYV + Y(k,k)*V(k)));
        end
        V(k) = (1/Y(k,k))*((P(k)-j*Q(k))/conj(V(k))-PYV); % Compute bus voltages
        if bus(k,10) == 2 % For PV buses, the voltage magnitude remains same, but the angle changes
            V(k)=abs(Vprev(k))*(cos(angle(V(k)))+j*sin(angle(V(k))));
        end
    end
    iteration=iteration+1; % Increment iteration count
    tolerance = max(abs(abs(V) - abs(Vprev))); % Tolerance at the current iteration
    Vprev=V;
end
ste=clock; % end the iteration time clock

%---------------------power flow----------------------%
% currents at each node
I=Y*V;
% Power at each node
S=V.*conj(I); % Complex power
for k=1:nbuses
    if bus(k,10)==1
        % Real and reactive generation at the Slack bus
        Pgen(k)=real(S(k));
        Qgen(k)=imag(S(k));
    end
    if bus(k,10)==2
        % Real and reactive generation at the PV buses
        Pgen(k)=real(S(k))+bus(k,6);
        Qgen(k)=imag(S(k))+bus(k,7);
    end
    if bus(k,10)==3
        Pgen(k)=0;
        Qgen(k)=0;
    end
end
% calculate the line flows and power losses
FromNode=line(:,1);
ToNode=line(:,2);
nbranch = length(line(:,1)); % number of branches
for k=1:nbranch
    a=line(k,6);   % Find out if is a line or a transformer, a=0 -> line, a=1 -> transformer, 0<a<1 -> Transformer
    switch a    % for both cases use the pi model
        case 0  %if its a line a=0
            b=1i*line(k,5);
            suceptancia(k,1)=b/2; 
            suceptancia(k,2)=b/2;
        otherwise %if its a transformer
            Zpq=line(k,3)+1i*line(k,4);
            Ypq=Zpq^-1;
            suceptancia(k,1)=(Ypq/a)*((1/a)-1); 
            suceptancia(k,2)=Ypq*(1-(1/a));
    end
end
% Define admmitance of lines
r = line(:,3);
rx = line(:,4);
z = r + j*rx;
y = ones(nbranch,1)./z;
% Define complex power flows
Ss = V(FromNode).*conj((V(FromNode) - V(ToNode)).*y ...
   + V(FromNode).*suceptancia(:,1)); % complex flow of the sending buses
Sr = V(ToNode).*conj((V(ToNode) - V(FromNode)).*y ...
   + V(ToNode).*suceptancia(:,2)); % complex low of the receiving buses

% Define active and reactive power flows
Pij=real(Ss);
Qij=imag(Ss);
Pji=real(Sr);
Qji=imag(Sr);

% Active power lossess
P_loss=sum(Pij+Pji);

% Reactive power lossess
Q_loss=sum(Qij+Qji);

% ----------------------Displaying the outputs----------%
disp('                      Gauss Seidel Load-Flow Study');
disp('                    Report of Power Flow Calculations ');
disp(' ');
disp(date);
fprintf('Number of iterations       : %g \n', iteration);
fprintf('Solution time              : %g sec.\n',etime(ste,st));
fprintf('Total real power losses    : %g.\n',P_loss);
fprintf('Total reactive power losses: %g.\n\n',Q_loss);
disp('                                      Generation             Load');
disp('      Bus      Volts     Angle      Real  Reactive      Real  Reactive ');
ywz=[  bus(:,1)    abs(V)  (180/pi)*angle(V)  Pgen'  Qgen'  bus(:,6)  bus(:,7)];
disp(ywz);
disp('                      Line Flows                     ');
disp('    #Line    From Bus   To Bus     Real    Reactive   ');
l=1:1:length(line(:,1));
xy=[l' FromNode ToNode Pij Qij];
yx=[l' ToNode  FromNode Pji Qji];
disp(xy);
disp(yx);

% --------------------------Y-bus building----------------%
 function [Y] = Y_admi(line,bus,nbuses)

orden=zeros(1,nbuses);
Y=zeros(nbuses);
for k=1:nbuses
    orden(k)=k; % vector which contains the order of building according to the bus information
end
for p=1:length(line(:,1));
    BusP=find(orden==line(p,1));
    BusQ=find(orden==line(p,2));
    a=line(p,6); %Tap value for the  p iteration
    if a>0 % for transformers out of nominal position
    yl=(1/((line(p,3)+j*line(p,4)))); % line admittance
    Ad=(j*line(p,5)/2); % line charging
    Y(BusP,BusQ)=Y(BusP,BusQ)-yl/a; % a non diagonal element
    Y(BusQ,BusP)=Y(BusP,BusQ); % symmetry is  declared for elements out of the diagonal
    Y(BusP,BusP)=Y(BusP,BusP)+(yl/a)+((1/a)*(1/a-1)*yl)+Ad; %Equivalent admittance at the P-terminal plus line charging
    Y(BusQ,BusQ)=Y(BusQ,BusQ)+(yl/a)+(1-1/a)*yl+Ad; %Equivalent admittance at the Q-terminal plus line charging
    else % for lines
         yl=(1/((line(p,3)+j*line(p,4)))); % line admittance
         Ad=(j*line(p,5)/2);  % line charging
         Y(BusP,BusQ)=Y(BusP,BusQ)-yl; % a non diagonal element
         Y(BusQ,BusP)=Y(BusP,BusQ); % symmetry is  declared for elements out of the diagonal
         Y(BusP,BusP)=Y(BusP,BusP)+yl; % diagonal element
         Y(BusQ,BusQ)=Y(BusQ,BusQ)+yl; % diagonal element
         c=line(p,5); % line charging for the whole line
         if c>0
             Y(BusP,BusP)= Y(BusP,BusP)+Ad; %add value of line charging to the diagonal element
             Y(BusQ,BusQ)= Y(BusQ,BusQ)+Ad; %add value of line charging to the diagonal element
         end

    end
end
% Susceptance and conductance
for m=1:nbuses
    dir=find(orden==bus(m,1));
    Y(dir,dir)=Y(dir,dir)+bus(m,8)+j*bus(m,9); % add the shunt admittance 
end
 end
 %-----------------------end-------------------------------%