clear all
close all 
syms d1 s q1 q2  l1 dq1 dq2 ddq1 ddq2  d1 d2 m1 m2 I1 I2 g real
L=[l1 l1];
q=[q1;q2];

%origins and centers
O0 = [0 0 0]';
Oc1 = [0 0 0]' ; %center of mass in B
O1 = [(l1)*cos(q1) (l1)*sin(q1) 0]';
Oc2 = [(l1+q2)*cos(q1) (l1+q2)*sin(q1) 0 ]' ; %q2 0 to L
O2 = [((2*q2)+l1)*cos(q1) ((2*q2)+l1)*sin(q1) 0 ]';

%Axes rotation
Z0 = [0 0 1]';
Z1 = [cos(q1) sin(q1) 0]';
Zer = [0 0 0]';


%Jacobbians
Jv1 = [cross(Z0,(Oc1-O0)), Zer];
Jv2 = [cross(Z0,(Oc2-O0)), Z1];
Jw1 = [Z0, Zer];
Jw2 = [Z0, Zer];
 

Vc1=Jv1(1:3,1)*dq1 + Jv1(1:3,2)*dq2;
Vc2=Jv2(1:3,1)*dq1 +  Jv2(1:3,2)*dq2;
Wc1=Jw1(1:3,1)*dq1 +  Jw1(1:3,2)*dq2;
Wc2=Jw2(1:3,1)*dq1 +  Jw2(1:3,2)*dq2;

K=0.5*m1*Vc1'*Vc1 + 0.5*m1*Vc2'*Vc2 + 0.5*Wc1'*I1*Wc1 + 0.5*Wc2'*I1*Wc2;
K=simplify(K)


R1 = [
cos(q1) -sin(q1) 0
sin(q1) cos(q1) 0
 0 0 1];
 
R2 = [
 cos(q1) -sin(q1) 0
 sin(q1) cos(q1) 0
 0 0 1];

 
%D
D1 = m1*Jv1'*Jv1 + Jw1'*R1*I1*R1'*Jw1;
D2 = m2*Jv2'*Jv2 + Jw2'*R2*I2*R2'*Jw2;
D=D1+D2;
D=simplify(D)
 

%P
P1 = m1*g*0*sin(q1);
P2 = m2*g*((l1+q2)*sin(q1));
P = P1+P2

%G
G1 = diff(P, q1);
G2 = diff(P, q2);
G = [G1; G2]

 
dq = [dq1; dq2];
ddq = [ddq1; ddq2];
%C
C = Coriolis(D,q,dq,2);
C=simplify(C)

%%
%%%%%%%%%%%%%%%%%%%%% Torques function
tor = D*ddq+C*dq+G;
tor=simplify(tor)
%%%%%%%%%%%%%%%%%%%%%
%%
%Appling a force
%subs
D(q1,q2) = subs(D,{m2, I1, I2, l1},{2 1 2 0.2});
C(q1,q2,dq1,dq2) = subs(C*dq ,{m2, I1, I2, l1},{2 1 2 0.2});
G(q1,q2) = subs(G ,{m2, I1, I2, l1,g},{2 1 2 0.2 9.8});
% 
%initial conditions
q1_0 = pi/2; % position 
q2_0 = 0;  
dq1_0 = 0; %velocity
dq2_0 = 0;
ddq1p_0 =0; %acceleration
ddq2p_0 =0;
dt=0.01;    %step in seconds
n=100;      %total steps

%force funtion to applied
u1p_0 = 10;
u2p_0 = 5;
for i = 1:n
   u1p(i)=u1p_0+(i*0.15);
   u2p(i)=u2p_0+(i*0.15); 
end

%Appling a force to the joints
for i = 1:n
     U = [u1p(i);u2p(i)];
     q1p(i)=q1_0;
     q2p(i)=q2_0;
     dq1p(i)=dq1_0;
     dq2p(i)=dq2_0;
     ddq1p(i)=ddq1p_0;
     ddq2p(i)=ddq2p_0;
     %inverse dinamics
     ddq = inv(D(q1_0, q2_0))*(U-C(q1_0, q2_0,dq1_0,dq2_0)-G(q1_0,q2_0));
     ddq1p_0 = ddq(1);
     ddq2p_0 = ddq(2);
     %small aceleration
     dq1_0=dq1p(i) + double(ddq(1)*dt);
     dq2_0=dq2p(i) + double(ddq(2)*dt);
     q1_0 = q1p(i) + dq1_0*dt;
     q2_0 = q2p(i) + dq2_0*dt;
end

%%
%Plots
t = 0:0.1:(0.1*(n-1));
 
 figure
 plot(t,u1p,'k','linewidth',2)
 hold on
 plot(t,u2p,'g','linewidth',2)
 grid on
 title('applied force(N) vs time(s)')
 legend('joint_1', 'joint_2')
 
 
 figure
 plot(t,q1p,'k','linewidth',2)
 hold on
 plot(t,q2p,'g','linewidth',2)
 grid on
 title('position(rad) vs time(s)')
 legend('joint_1', 'joint_2')
 
 figure
 plot(t,dq1p,'k','linewidth',2)
 hold on
 plot(t,dq2p,'g','linewidth',2)
 grid on
 title('velocity(rad/s) vs time(s)')
 legend('joint_1', 'joint_2')
 

 figure
 plot(t,ddq1p,'k','linewidth',2)
 hold on
 plot(t,ddq2p,'g','linewidth',2)
 grid on
 title('acceleration(rad/s^2) vs time(s)')
 legend('joint_1', 'joint_2')
 

 

%%
%functions
function C = Coriolis(D,q,dq,n)
    sym C;
    for k = 1:n
        for j =1:n
        C(k,j) = sym(0);
            for i=1:n
                c_ijk = (1/2)*(diff(D(k,j),q(i)) + diff(D(k,i),q(j))-diff(D(i,j),q(k)));
                C(k,j) = C(k,j) + c_ijk*dq(i);
            end
        end
    end
end