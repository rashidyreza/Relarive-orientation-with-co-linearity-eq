kkkk=input('k=')
switch kkkk
    case 1 
xp=0;yp=0;
[n1,n2]=uigetfile('*.txt','xyL');%xy photo grametric aks chap
XY_L=strcat(n2,n1);
XY_L=dlmread(XY_L);
XY_L=0.001.*XY_L;

[n1,n2]=uigetfile('*.txt','xyR');%xy photo grametric aks rast
XY_R=strcat(n2,n1);
XY_R=dlmread(XY_R);
XY_R=0.001.*XY_R;
n=size(XY_R,1);
Qy=eye(4*n);


%focal lens
f = .152844;
syms omega2 phi2 kapa2 by bz %Unknown
%Rotation matrix
R_kapa = [cos(kapa2) sin(kapa2) 0 ; -sin(kapa2) cos(kapa2) 0 ; 0 0 1];
R_phi = [cos(phi2) 0 -sin(phi2) ; 0 1 0 ; sin(phi2) 0 cos(phi2)];
R_omega = [1 0 0 ; 0 cos(omega2) sin(omega2) ; 0 -sin(omega2) cos(omega2)];
R = R_kapa*R_phi*R_omega;

xL = sym('xL',[1 n]);
yL = sym('yL',[1 n]);
xR = sym('xR',[1 n]);
yR = sym('yR',[1 n]);

bx = 1;
% Coplanarity Equation
for i=1:n
 F0(i,1)=det( [bx by bz ; xL(i)-xp yL(i)-yp -f ;[xR(i)-xp yR(i)-yp -f]*R'] );
end
L = [xL yL xR yR];
L0=[XY_L(:,1)',XY_L(:,2)',XY_R(:,1)',XY_R(:,2)'];
L00=L0;
unknown = [omega2 phi2 kapa2 by bz];
B0 = jacobian(F0,L);
A0 = jacobian(F0,unknown);
unknown0 = [0 0 0 0 0];
e0=2;
e1=3;
dx=1;
while (((norm(dx)< 1e-8)&&(norm(e1-e0))<1e-8))==0
    
    A = subs(A0,L,L0);
    F = subs(F0,L,L0);
    B = subs(B0,L,L0);
    A = eval(subs(A,unknown,unknown0));
    B = eval(subs(B,unknown,unknown0));
    W = eval(subs(F,unknown,unknown0))+B*(L00'-L0');
    dx = -inv(A'*inv(B*Qy*B')*A)*A'*inv(B*Qy*B')*W;
    e0=e1;
    e1=Qy*B'*inv(B*Qy*B')*(eye(n)-A*inv(A'*inv(B*Qy*B')*A)*A'*inv(B*Qy*B'))*W;
    
    unknown0 = unknown0 + dx';
    L0=L00-e1';
end
omega2=unknown0(1);phi2=unknown0(2);kapa2=unknown0(3);by=unknown0(4);bz=unknown0(5);
disp('     omega2                phi2                    kapa2                  by                 bz')
unknown0


 xxx=[XY_L(:,1)-xp XY_L(:,2)-yp -f*ones(n,1)];
x_1=xxx(:,1);
y_1=xxx(:,2);
z_1=xxx(:,3);
 XXX=[XY_R(:,1)-xp XY_R(:,2)-yp -f*ones(n,1)]*eval(R');
x_2=XXX(:,1);
y_2=XXX(:,2);
z_2=XXX(:,3);
landa=(bx.*z_2-bz.*x_2)./(x_1.*z_2-x_2.*z_1);
mu=(bx.*z_1-bz.*x_1)./(x_1.*z_2-x_2.*z_1);
%Model Coordinates
Xm=landa.*x_1
Zm=landa.*z_1
Ym=0.5*landa.*y_1+0.5*mu.*y_2+0.5*by
py=landa.*y_1-mu.*y_2-by
dlmwrite('XYZ_Model.txt',[Xm,Ym,Zm],'newline','pc','precision',12);

    case 2 
[n1,n2]=uigetfile('*.txt','xyL');%xy photo grametric aks chap
XY_L=strcat(n2,n1);
XY_L=dlmread(XY_L);

[n1,n2]=uigetfile('*.txt','xyR');%xy photo grametric aks rast
XY_R=strcat(n2,n1);
XY_R=dlmread(XY_R);
n=size(XY_R,1);
Qy=eye(4*n);

%focal lens

f = 152.844;

syms omega2 phi2 kapa2 kapa1 phi1 %Unknown

%Rotation matrix
R2_kapa = [cos(kapa2) sin(kapa2) 0 ; -sin(kapa2) cos(kapa2) 0 ; 0 0 1];
R2_phi = [cos(phi2) 0 -sin(phi2) ; 0 1 0 ; sin(phi2) 0 cos(phi2)];
R2_omega = [1 0 0 ; 0 cos(omega2) sin(omega2) ; 0 -sin(omega2) cos(omega2)];
R2 = [R2_kapa*R2_phi*R2_omega]';
R1_kapa = [cos(kapa1) sin(kapa1) 0 ; -sin(kapa1) cos(kapa1) 0 ; 0 0 1];
R1_phi = [cos(phi1) 0 -sin(phi1) ; 0 1 0 ; sin(phi1) 0 cos(phi1)];
R1 = [R1_kapa*R1_phi]';

xL = sym('xL',[1 n]);
yL = sym('yL',[1 n]);
xR = sym('xR',[1 n]);
yR = sym('yR',[1 n]);

bx = 1;
% Coplanarity Equation
for i=1:n
 F0(i,1)=det( [bx 0 0 ; [xL(i) yL(i) -f]*R1' ;[xR(i) yR(i) -f]*R2'] );
end
L = [xL yL xR yR];
L0=[XY_L(:,1)',XY_L(:,2)',XY_R(:,1)',XY_R(:,2)'];
L00=[XY_L(:,1)',XY_L(:,2)',XY_R(:,1)',XY_R(:,2)'];
unknown = [omega2 phi2 kapa2 kapa1 phi1];
B0 = jacobian(F0,L);
A0 = jacobian(F0,unknown);
unknown0 = [0 0 0 0 0];


dx=1;
e0=2;
e1=3;
while (((norm(dx)< 1e-12)&&(norm(e1-e0))<1e-12))==0
    A = subs(A0,L,L0);
    F = subs(F0,L,L0);
    B = subs(B0,L,L0);
    A = eval(subs(A,unknown,unknown0));
    B = eval(subs(B,unknown,unknown0));
    W = eval(subs(F,unknown,unknown0))+B*(L00'-L0');
    dx = -inv(A'*inv(B*Qy*B')*A)*A'*inv(B*Qy*B')*W;
    e0=e1;
    e1=Qy*B'*inv(B*Qy*B')*(eye(n)-A*inv(A'*inv(B*Qy*B')*A)*A'*inv(B*Qy*B'))*W;
    unknown0 = unknown0 + dx';
    L0=L00-e1';
end
disp('     omega2                phi2                    kapa2                  kapa1                 phi1')
omega2=unknown0(1);phi2=unknown0(2);kapa2=unknown0(3);kapa1=unknown0(4);phi1=unknown0(5);
unknown0
end














% rL = [xl; yl; zl];
% rR = [xr; yr; zr];
% 
% bb = [0 bz -by;-bz 0 bx;by -bx 0];
% 
% MWL = [1 0 0;0 cos(WL) -sin(WL);0 sin(WL) cos(WL)];
% MPL = [cos(PL) 0 sin(PL);0 1 0;-sin(PL) 0 cos(PL)];
% MKL = [cos(KL) -sin(KL) 0;sin(KL) cos(KL) 0;0 0 1];
% 
% MWR = [1 0 0;0 cos(WR) -sin(WR);0 sin(WR) cos(WR)];
% MPR = [cos(PR) 0 sin(PR);0 1 0;-sin(PR) 0 cos(PR)];
% MKR = [cos(KR) -sin(KR) 0;sin(KR) cos(KR) 0;0 0 1];
% 
% MML = MKL*MPL*MWL;
% MMR = MKR*MPR*MWR;
% 
% ff = transpose(rL)*MML*bb*transpose(MMR)*rR;