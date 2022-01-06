clear all;
close all;
clc;

%-----------------All values are in SI units------------------------%
global g rho_f vis_f h u0 d_p rho_p dt
% Gravity
g=-9.81;
% Fluid property(Air)
rho_f = 1.2;
vis_f = 1.8 * 10^-5;
u0 =0.2 ;
% fluid flow height 
h = 0.1;
rho_g = 750;      
d_g = 2.5 * 10^-3;
rho_c = 50;       
d_c =3.25 * 10^-3;
x0 = 0.5;  y0 = 0.5; up0= 0; vp0 = 0;
x_c =0.55; % bin separation point 
y_c = -0.5;
%-------------------------case [1]---------------------------------%

 %-----------------------------Euler method----------------------%

  %------------------------------grain--------------------------%
rho_p=rho_g ;
d_p =d_g ;
x=x0;
y=y0;
u=up0;
v=vp0;
dt=1e-2;
ii=1;
count = 0;
while(y(ii)>=y_c)
    [du1,dv1] = sys_01(x(ii),y(ii),u(ii),v(ii));
    x(ii+1)=x(ii)+ dt*(u(ii));
    y(ii+1)=y(ii)+ dt*(v(ii));
    u(ii+1) = u(ii)+dt*(du1);
    v(ii+1)= v(ii)+dt*(dv1);
    ii =ii+1;
end
X_euler_gr = x(end);

count;
figure(1)
plot(x,y,'g')
hold on;

x_euler_g = x;
y_euler_g = y;
u_euler_g = u;
v_euler_g = v;

clear x y u v

%  -------------------------chaff------------------------------%

rho_p=rho_c ;
d_p =d_c ;
x=x0;
y=y0;
u=up0;
v=vp0;
y_c = -0.5;
dt=1e-2;
Flag =0;
ii =1;
while (y(ii)>=y_c)
    [du1,dv1] = sys_01(x(ii),y(ii),u(ii),v(ii));
    x(ii+1)=x(ii)+ dt*(u(ii));
    y(ii+1)=y(ii)+ dt*(v(ii));
    u(ii+1) = u(ii)+dt*(du1);
    v(ii+1)= v(ii)+dt*(dv1);
    ii= ii+1;
    Flag =Flag+1;
end
X_euler_ch = x(end);

Flag;
figure(1)
plot(x,y,'r')
hold off;
title('Euler Method');
xlabel('X-distance')
ylabel('Y-diatance')
legend('grain particle','chaff particle')

x_euler_c = x;
y_euler_c = y;
u_euler_c = u;
v_euler_c = v;


clear x y u v

%----------------------------Runga-kuta-4---------------------%

% grains
rho_p=rho_g ;
d_p =d_g ;
x=x0;
y=y0;
u=up0;
v=vp0;
y_c = -0.5;
dt=1e-2;
ii=1;
rounds =0;
while (y(ii)>=y_c)
    [x(ii+1),y(ii+1),u(ii+1),v(ii+1)]=sys_02(x(ii),y(ii),u(ii),v(ii));
    ii= ii+1;
end 
X_Rk4_gr = x(end);
rounds = rounds +1;
figure(2)
plot(x,y,'b')
hold on;
x_Rk_g = x;
y_Rk_g = y;
u_Rk_g = u;
v_Rk_g = v;

clear x y u v

% chaff
rho_p=rho_c ;
d_p =d_c ;
x=x0;
y=y0;
u=up0;
v=vp0;
y_c = -0.5;
dt=1e-2;
ii=1;
rounds =0;
while (y(ii)>=y_c)
    [x(ii+1),y(ii+1),u(ii+1),v(ii+1)]=sys_02(x(ii),y(ii),u(ii),v(ii));
    ii= ii+1;
end   
X_Rk4_ch = x(end);

rounds = rounds +1;
figure(2)
plot(x,y,'g')
hold off 

title('RK4 Method');
xlabel('X-distance')
ylabel('Y-diatance')
legend('grain particle','chaff particle')

x_Rk_c = x;
y_Rk_c = y;
u_Rk_c = u;
v_Rk_c = v;

clear x y u v



%---------------------------case[2]--------------------------%
%----------------------euler accuracy------------------------%

% --------------grain and chaff particle --------------------%
%grain
rho_p=rho_g ;
d_p =d_g ;
x=x0;
y=y0;
u=up0;
v=vp0;
y_c = -0.5;
dt=1e-5;
ii=1;
count = 0;
while(y(ii)>=y_c)
    [du1,dv1] = sys_01(x(ii),y(ii),u(ii),v(ii));
    x(ii+1)=x(ii)+ dt*(u(ii));
    y(ii+1)=y(ii)+ dt*(v(ii));
    u(ii+1) = u(ii)+dt*(du1);
    v(ii+1)= v(ii)+dt*(dv1);
    ii =ii+1;
end
X_euler_accuracy_gr = x(end);

figure(3)
subplot(2,2,1)
plot(x,y,'b-',x_euler_g,y_euler_g,'r.-')
title('Euler accuracy grain');
xlabel('X-distance')
ylabel('Y-diatance')
legend('grain particle','chaff particle')


%chaff
rho_p=rho_c ;
d_p =d_c ;
x=x0;
y=y0;
u=up0;
v=vp0;
y_c = -0.5;
dt=1e-5;
Flag =0;
ii =1;
while (y(ii)>=y_c)
    [du1,dv1] = sys_01(x(ii),y(ii),u(ii),v(ii));
    x(ii+1)=x(ii)+ dt*(u(ii));
    y(ii+1)=y(ii)+ dt*(v(ii));
    u(ii+1) = u(ii)+dt*(du1);
    v(ii+1)= v(ii)+dt*(dv1);
    ii= ii+1;
    Flag =Flag+1;
end
X_euler_accuracy_ch = x(end);

figure(3)
subplot(2,2,2)
plot(x,y,'k-',x_euler_c,y_euler_c,'g.-')
title('Euler accuracy chaff ');
xlabel('X-distance')
ylabel('Y-diatance')
legend('grain particle','chaff particle')


%----------------------Rk4 method accuracy------------------------%

   % --------------grain and chaff particle --------------------%
% grains
rho_p=rho_g ;
d_p =d_g ;
x=x0;
y=y0;
u=up0;
v=vp0;
y_c = -0.5;
dt=1e-4;
ii=1;
rounds =0;
while (y(ii)>=y_c)
    [x(ii+1),y(ii+1),u(ii+1),v(ii+1)]=sys_02(x(ii),y(ii),u(ii),v(ii));
    ii= ii+1;
end 
X_Rk4_accuracy_gr = x(end);

figure(3)
subplot(2,2,3)
plot(x,y,'k-',x_Rk_g,y_Rk_g,'g.-')
title('RK4 accuracy grain');
xlabel('X-distance')
ylabel('Y-diatance')
legend('grain particle','chaff particle')

% chaff
rho_p=rho_c ;
d_p = d_c ;
x=x0;
y=y0;
u=up0;
v=vp0;
y_c = -0.5;
dt=1e-4;
ii=1;
rounds =0;
while (y(ii)>=y_c)
    [x(ii+1),y(ii+1),u(ii+1),v(ii+1)]=sys_02(x(ii),y(ii),u(ii),v(ii));
    ii= ii+1;
end
X_Rk4_accuracy_ch = x(end);

figure(3)
subplot(2,2,4)
plot(x,y,'b-',x_Rk_c,y_Rk_c,'r.-')
title('RK4 accuracy chaff');
xlabel('X-distance')
ylabel('Y-diatance')
legend('grain particle','chaff particle')
set(gcf,'position',[10 10 2000 2000])


% accuracy 
% for euler method (grain)
D_error_euler_g = abs(X_euler_accuracy_gr - X_euler_gr )/X_euler_accuracy_gr
D_error_euler_c = abs(X_euler_accuracy_ch - X_euler_ch)/X_euler_accuracy_ch
D_error_Rk4_g = abs(X_Rk4_accuracy_gr - X_Rk4_gr)/X_Rk4_accuracy_gr
D_error_Rk4_c = abs(X_Rk4_accuracy_ch - X_Rk4_ch)/X_Rk4_accuracy_ch



%---------------------------case[3]------------------------%

%  varing jet velocity to make sure , chaff particle fall inside the bin 2
n = 1;
rho_c = 50;       
d_c =3.25 * 10^-3; % chaff particle property 
dt =1e-2;
while true
    u0 = u0-0.001;
    clear x y u v 
    x = x0;
    y = y0;
    u = up0;
    v = up0;
    ii=1;
    while (y(ii)>=y_c)
        [x(ii+1),y(ii+1),u(ii+1),v(ii+1)]=sys_02(x(ii),y(ii),u(ii),v(ii));
        ii= ii+1;
    end 
    if (x(end)<x_c)
        break;
    end 
end 
min_velocity = u0+0.001

    %-----------------------case[4]------------------------%
    
 % Derive terminal velocity 
 Ut_g = sqrt(d_g*(rho_f-rho_g)*g/(3*rho_f*0.44)); %grain (heavy)
 Ut_c = sqrt(d_c*(rho_f-rho_c)*g/(3*rho_f*0.44)); %chaff (light)
 angle_1 = deg2rad(-108);
 angle_2 = deg2rad(-72);
 u0 = 0.2;
 
 
% part_a (change x_c)

Utr_g_x = Ut_g*cos(angle_2);
Utr_C_x = Ut_c*cos(angle_1);
Utr_g_y = Ut_g*sin(angle_2);
Utr_c_y = Ut_c*sin(angle_1);
up = [Utr_g_x, Utr_C_x];
vp = [Utr_g_y, Utr_c_y];
y_c = -0.5;
dt=1e-2;

rho=[rho_g,rho_c];
d = [d_g,d_c];
for j=1:length(rho)
    ii=1;
    clear x y u v 
    x=x0;
    y=y0;
    u=up(j);
    v=vp(j);
    rho_p=rho(j) ;
    d_p =d(j) ;
    while (y(ii)>=y_c)
        [x(ii+1),y(ii+1),u(ii+1),v(ii+1)]=sys_02(x(ii),y(ii),u(ii),v(ii));
        ii= ii+1;
    end 
    tra_x(j)=x(end);
end

if (tra_x(1) < tra_x(2)) % (grain trajectory < chaff trajectory)
    disp('We can predict a complete separation');
else 
    disp('There are chaff particle in bin 1');
end 

% part c (change u0)

Utr_g_x = Ut_g*cos(angle_2);
Utr_C_x = Ut_c*cos(angle_1);
Utr_g_y = Ut_g*sin(angle_2);
Utr_c_y = Ut_c*sin(angle_1);
up = [Utr_g_x, Utr_C_x];
vp = [Utr_g_y, Utr_c_y];
x_c = 0.55;
y_c =-0.5;
dt=1e-2;
rho=[rho_g,rho_c];
d = [d_g,d_c];

for u0 = 0.05:0.001:0.4
    for j = 1:length(d)
    ii=1;
    clear x y u v 
    x=x0;
    y=y0;
    u=up(j);
    v=vp(j);
    rho_p=rho(j) ;
    d_p =d(j) ;
    while (y(ii)>=y_c)
        [x(ii+1),y(ii+1),u(ii+1),v(ii+1)]=sys_02(x(ii),y(ii),u(ii),v(ii));
        ii= ii+1;
    end 
    tra_x(j)=x(end);
    end 
end 

if (tra_x(1) < tra_x(2)) % (grain trajectory < chaff trajectory)
    disp('We can predict a complete separation');
else 
    disp('There are chaff particle in bin 1');
end 



% --------------defining functions--------------%

function[x2,y2,u2,v2] = sys_02(x1,y1,u1,v1)
global dt 
 [du1,dv1] = sys_01(x1,y1,u1,v1);
    k1_x = u1;
    k1_y = v1;
    k1_u = du1;
    k1_v = dv1;
    
    [du2,dv2] = sys_01(x1+dt*u1/2,y1+ dt*v1/2,u1+dt*du1/2,v1+dt*dv1/2);
    k2_x = u1+ du1*dt /2;
    k2_y = v1+ dv1*dt /2;
    k2_u = du2;
    k2_v = dv2;
    
    [du3,dv3] = sys_01(x1+dt*u1/2,y1+ dt*v1/2,u1+dt*du2/2,v1+dt*dv2/2);
    k3_x = u1+ du2*dt /2;
    k3_y = v1+ du2*dt /2;
    k3_u = du3;
    k3_v = dv3;
    
    [du4,dv4] = sys_01(x1+dt*u1,y1+ dt*v1,u1+dt*du3,v1+dt*dv3);
    k4_x = u1+ du3*dt;
    k4_y = v1+ du3*dt;
    k4_u = du4;
    k4_v = dv4;
    
    x2 = x1+ dt*(1/6)*(k1_x+2*k2_x+2*k3_x+k4_x);
    y2 = y1+ dt*(1/6)*(k1_y+2*k2_y+2*k3_y+k4_y);
    u2 = u1+ dt*(1/6)*(k1_u+2*k2_u+2*k3_u+k4_u);
    v2 = v1+ dt*(1/6)*(k1_v+2*k2_v+2*k3_v+k4_v);
end 
% function for derivative of velocity.

function [F_tu,F_tv]= sys_01(x,y,up,vp)  
global g rho_f vis_f h d_p rho_p u0
% Equations
uf = 6.2*u0*exp(-50*((y/x)^2))*(h/x)^(1/2);
vf = 0;
V_p = (4/3)*pi*(d_p/2)^3;
F_gravity = rho_p*V_p*g;
F_buoyancy = -rho_f*V_p*g;

% since F_drag_x and y depends Re_p and Cd
Re_p = rho_f* sqrt((uf-up)^2 + (vf-vp)^2) * d_p/vis_f;
if Re_p < 800
    Cd = (24/Re_p)*(1+(0.15*(Re_p)^(0.687)));
else 
    Cd = 0.44;
end 
F_drag_x = (0.5*pi*rho_f*(d_p^2)*Cd*sqrt((uf-up)^2 + (vf-vp)^2)*(uf-up));
F_drag_y = (0.5*pi*rho_f*(d_p^2)*Cd* sqrt((uf-up)^2 + (vf-vp)^2)*(vf-vp));

F_tu = F_drag_x/(rho_p*V_p);
F_tv = ( F_gravity+F_buoyancy+F_drag_y )/(rho_p*V_p);
end 
   