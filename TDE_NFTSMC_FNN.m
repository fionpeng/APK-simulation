function [sys,x0,str,ts] = tdesmc(t,x,u,flag)

switch flag,
case 0,
    [sys,x0,str,ts]=mdlInitializeSizes;
case 2,
    sys=mdlUpdate(t,x,u);
case 3,
    sys=mdlOutputs(t,x,u);
case {4,9}
    sys=[];
otherwise
    error(['Unhandled flag = ',num2str(flag)]);
end

function [sys,x0,str,ts]=mdlInitializeSizes
sizes = simsizes;
sizes.NumOutputs     = 6;
sizes.NumInputs      = 2;
sizes.NumDiscStates = 134;
sizes.DirFeedthrough = 1;
sizes.NumSampleTimes = 1;
sys = simsizes(sizes);
x0  = zeros(1,134);
str = [];
ts  = [0.001 0];

function sys=mdlUpdate(t,x,u) %
persistent ka
if t==0
    ka=0;
end
a=5
M=0.001;
L=0.001;
k=1;   %
fn=1;
lamba=0.5;
%FNN

alfa=0.05;
xite=0.35;
%FNN

bi=[2 2 2 2 2];
ci=[-3 -1 0 1 3;-3 -1 0 1 3;-3 -1 0 1 3];

qd=u(1);
dqd=15*pi*cos(pi*t-pi/2);
ddqd=-15*pi*pi*sin(pi*t-pi/2);
p=3;
l=5;
fi=2;
the1=1; %调节后误差变化明显
the2=0.1; %调节后误差变化明显
enta=6;
T=4;
ex=0.01;

q=u(2);
dq=(q-x(1))/L;
ddq=(dq-x(2))/L;
% ddq=(u(2)-2*q_xx+q_x)/L^2;

e=q-qd;
de=dq-dqd;
dde=ddq-ddqd
ee=x(133)+e;
%M.JIN2015 no used for comparing
aaf=1;
bbf=1;
kkf=0.1;
f1=0.1;
f2=1;

% bei=0.5;
% s=e+the1*e^fi+the2*de^(l/p)
% s=e+bei*abs(de)^(l/p)*sign(de); %FTSMC M.JIN2015  for comparing

s=e+the1*abs(e)^fi*sign(e)+the2*abs(de)^(l/p)*sign(de) %NFTSM sliding surface
H=x(4)-M*x(3);
;
%FNN
xi=[e,de,dde]';
FS1=0;
for j=1:1:5
   u1(j)=exp(-norm(xi(1)-ci(1,j))^2/(2*bi(j)*bi(j)));
end
for j=1:1:5
     gs2=-norm(xi(2)-ci(2,j))^2/(2*bi(j)*bi(j));
     u2(j)=exp(gs2);
end
for j=1:1:5
    u3(j)=exp(-norm(xi(3)-ci(3,j))^2/(2*bi(j)*bi(j)));
end
% 

for l1=1:1:5
	for l2=1:1:5
        for l3=1:1:5
            FS2(5*5*(l1-1)+5*(l2-1)+l3)=u1(l1)*u2(l2)*u3(l3);
            FS1=FS1+u1(l1)*u2(l2)*u3(l3);
        end
	end
end
FS=FS2/(FS1+0.001);
kk=0.1; %
gama=2;
SS=(1/M)*s*the2*(l/p)*abs(de)^(l/p-1)*FS*kk;  %FNN Formula (28)
for i=1:1:5*5*5
    thta(i,1)=x(i+4);
end
unn=-thta'*FS'


% KD=0.05;KP=0.02   
% s=de+KD*e+KP*ee;
% tu=H+M*(ddqd+KD*de+KP*e); %conventional formulation of the TDC iPD controller
% KP=0.05;KI=0.0001;
% tu=H+M*(ddqd+KI*ee+KP*e); %iPI of the TDC
% KP=0.02;KI=0.01;KD=0.05;
% tu=H+M*(ddqd+KP*e+KI*ee+KD*de);
% tu=H+M*(ddqd+k*de+slaw);  %指数趋近滑模控制
% ueq=H+M*(ddq+lamba*(l/p)*abs(de)^(l/p-1)*de)
T1=abs(de)^(2-l/p)*sign(de);
ueq=-M*1/the2*p/l*(T1+the1*fi*abs(e)^(fi-1)*T1)+H+M*ddqd %NFTSMC
% ueq=-M*1/bei*p/l*T1+H+M*ddqd % FTSMC   M.JIN2015 
% ueq= -M*bei*l/p*abs(e)^(l/p-1)*sign(e)*de+H+M*ddqd  %FTSMC,error

if s< -ex
    sat=-1;
elseif s>ex
        sat=1;
else sat=s/ex;
end


bb=20;
ka=ka+x(133)*L
 ucor=-M*(ka+bb)*sat; %有自适应率
% ucor=-M*(10*s+1*sign(s)*abs(s)^(p/l));  %M.JIN2015
% ucor=-M*enta/T*sat; %where η> 0 and T> 0 are the convergence factor and the boundary layer thickness, respectively
dka=s*the2*(l/p)*fi*abs(de)^(l/p-1)*sign(s)  %鲁棒性比普通的要好 
tu=ueq+unn+ucor%NFTSM

% s=e+lamba*abs(de)^(l/p)*sign(de) ;%FNSMC M.JIN2015
% tu=M*(ddq+1/(lamba*l/p)*abs(de)^(2-l/p)*sign(de))-M*(the1*s+the2*abs(s)^(p/l)*sign(s))+H;
% tu=ueq+ucor

%  if t==2                      % Td
%     tu=tu+7;
% end 


% td1=8*q+10*dq^2+5*cos(q)       % Td1 and Td2
% if t==4
%     tu=td1+tu;
% end 
% 
% if t>1.5&&t<1.53
%     tu=0.5*tu;
% end 

sys(1)=u(2);
sys(2)=dq;
sys(3)=ddq;
sys(4)=tu;
for i=5:1:129
    sys(i)=SS(i-4);
end
sys(130)=dqd;
sys(131)=ddqd;
sys(132)=s;
sys(133)=dka;
sys(134)=H;


function sys=mdlOutputs(t,x,u)

sys(1)=x(4); %u
sys(2)=x(2); %dq
sys(3)=x(132); %s
sys(4)=x(130); %dqd
sys(5)=x(131); %ddqd
sys(6)=x(134); %H
