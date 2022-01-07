%--------------------------Assignment 2---------------------------------%
%------------------Srijith Sreekumar "230154" ----------------------------%
    close all;
    clc;
    tic

% set random number generator to default 
rng(5);% rng(seed) seeds the random number generator using the nonnegative integer seed 
%so that rand, randi, and randn produce a predictable sequence of numbers.


%%%%%%%% All the function are mentioned at end of the script file%%%%%%%%

%-----------------All values are in SI units------------------------%
global g rho_f vis_f h u0 d_p rho_p dt N

% Gravity
g=-9.81;%[m/s^2]

% Fluid property(Air)
rho_f = 1.2;%[kg/m^3]
vis_f = 1.8 * 10^-5;%[Pas]
u0 =0.2 ;%[m/s]

% fluid flow height 
h = 0.1;%[m]

% particle data
rho_g = 750;%[kg/m^3]      
d_g = 2.5 * 10^-3;%[m]
rho_c = 50; %[kg/m^3]      
d_c =3.25 * 10^-3;%[m]
x0 = 0.5;  y0 = 0.5; up0= 0; vp0 = 0;%[m,m.m/s,m/s] respectively
x_c =0.55; % [m]bin separation point 
y_c = -0.5;%[m]

%-------------------------------case [1]---------------------------%

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
while(y(ii)>=y_c)
    [du1,dv1] = particle_velocity_ME(x(ii),y(ii),u(ii),v(ii));
    x(ii+1)=x(ii)+ dt*(u(ii));
    y(ii+1)=y(ii)+ dt*(v(ii));
    u(ii+1) = u(ii)+dt*(du1);
    v(ii+1)= v(ii)+dt*(dv1);
    ii =ii+1;
end

figure(1)
plot(x,y,'g')
hold on;

x_euler_g = x;
y_euler_g = y; % the valuses of the x,y,u,v is assigned to other variables for further usage 
u_euler_g = u;
v_euler_g = v;

clear x y u v % removes the stored values in these varaibles.

%  -------------------------chaff------------------------------%

rho_p=rho_c ;
d_p =d_c ;
x=x0;
y=y0;
u=up0;
v=vp0;
y_c = -0.5;
dt=1e-2;
ii =1;
while (y(ii)>=y_c)
    [du1,dv1] = particle_velocity_ME(x(ii),y(ii),u(ii),v(ii));
    x(ii+1)=x(ii)+ dt*(u(ii));
    y(ii+1)=y(ii)+ dt*(v(ii));
    u(ii+1) = u(ii)+dt*(du1);
    v(ii+1)= v(ii)+dt*(dv1);
    ii= ii+1;
end

figure(1)
plot(x,y,'r')
plot(x_c,y_c,'*k')
hold off;
title('Euler Method');
xlabel('X-distance [m]')
ylabel('Y-diatance [m]')
xlim([0.49,0.64])
legend('grain particle[Euler]','chaff particle[Euler]','separation point')
saveas(gcf,'Euler profile[G&C].jpg')

x_euler_c = x;
y_euler_c = y;
u_euler_c = u;
v_euler_c = v;

clear x y u v 

%-------------------------Monte-carlo simulation-----------------------%
    %------------------------------grain--------------------------%
rho_p=rho_g ;
d_p =d_g ;
N = 100;
x=x0;
y=y0;
u=up0;
v=vp0;
dt=1e-2;
ii=1;
while(y(ii)>=y_c)
    % The monte carlo appraoach is written as a finction 
    [x(ii+1),y(ii+1),u(ii+1),v(ii+1)] = monte_carlo(x(ii),y(ii),u(ii),v(ii));
    ii = ii+1;
end 
  
figure(2)
plot(x,y,'g')
hold on;

x_monte_g = x;
y_monte_g = y;
u_monte_g = u;
v_monte_g = v;

clear x y u v

    %-------------------------chaff------------------------------%
rho_p=rho_c ;
d_p =d_c ;
N =100;
x=x0;
y=y0;
u=up0;
v=vp0;
y_c = -0.5;
dt=1e-2;
ii =1;

while(y(ii)>=y_c)
    [x(ii+1),y(ii+1),u(ii+1),v(ii+1)] = monte_carlo(x(ii),y(ii),u(ii),v(ii));
    ii = ii+1;
end 

figure(2)
plot(x,y,'r')
plot(x_c,y_c,'*k')
hold off;
title('Monte-Carlo Method');
xlabel('X-distance [m]')
ylabel('Y-diatance [m]')
xlim([0.49,0.64])
legend('grain particle[MC]','chaff particle[MC]','separation point')
saveas(gcf,'Monte-carlo profile[G&C].jpg')

x_monte_c = x;
y_monte_c = y;
u_monte_c = u;
v_monte_c = v;

clear x y u v

%-------------------------case 2---------------------------------------%
 %-----------------------Accuracy appraoch----------------------------%
 %    [ The step dt is further reduced to check the accuracy. ]
 
 %-----------------------Euler accuracy [grains]----------------------%
rho_p=rho_g ;
d_p =d_g ;
x=x0;
y=y0;
u=up0;
v=vp0;
dt=1e-4;
ii=1;
while(y(ii)>=y_c)
    [du1,dv1] = particle_velocity_ME(x(ii),y(ii),u(ii),v(ii));
    x(ii+1)=x(ii)+ dt*(u(ii));
    y(ii+1)=y(ii)+ dt*(v(ii));
    u(ii+1) = u(ii)+ dt*(du1);
    v(ii+1)= v(ii)+ dt*(dv1);
    ii =ii+1;
end
X_euler_gr = x(end);

x_euler_accuracy_g = x;
y_euler_accuracy_g = y;  
u_euler_accuracy_g = u;
v_euler_accuracy_g = v;
 
figure(3)
subplot (2,2,1)
plot(x_euler_g,y_euler_g,'g')
hold on;
plot(x_euler_accuracy_g,y_euler_accuracy_g,'k-')
%plot(x_c,y_c,'*k')
hold off
title("Euler accuracy [Grain]");
xlabel("X distance [m]");
ylabel("Y distance [m]");
%xlim([0.49,0.56])
legend("Grain Particles","Accuracy of Grain Particles");

clear x y u v

%  -----------------Euler accuracy [chaff]--------------------------%

rho_p=rho_c ;
d_p =d_c ;
x=x0;
y=y0;
u=up0;
v=vp0;
y_c = -0.5;
dt=1e-4;
ii =1;
while (y(ii)>=y_c)
    [du1,dv1] = particle_velocity_ME(x(ii),y(ii),u(ii),v(ii));
    x(ii+1)=x(ii)+ dt*(u(ii));
    y(ii+1)=y(ii)+ dt*(v(ii));
    u(ii+1) = u(ii)+dt*(du1);
    v(ii+1)= v(ii)+dt*(dv1);
    ii= ii+1;
end
X_euler_ch = x(end);

x_euler_accuracy_c = x;
y_euler_accuracy_c = y;
u_euler_accuracy_c = u;
v_euler_accuracy_c = v;

figure(3)
subplot(2,2,2)
plot(x_monte_c,y_monte_c,'r')
hold on
plot(x_euler_accuracy_c,y_euler_accuracy_c,'k-')
plot(x_c,y_c,'*k')
hold off
title("Euler Accuracy [Chaff] ");
xlabel("X distance [m]");
ylabel("Y distance [m]");
xlim([0.49,0.64])
legend("Chaff Particles","Accuracy of Chaff Particles",'separation point');

clear x y u v 

%-------------------Monte-carlo accuracy [grains]-----------------------%
    
rho_p=rho_g ;
d_p =d_g ;
N =100;
x=x0;
y=y0;
u=up0;
v=vp0;
dt=1e-4;
ii=1;
while(y(ii)>=y_c)
    [x(ii+1),y(ii+1),u(ii+1),v(ii+1)] = monte_carlo(x(ii),y(ii),u(ii),v(ii));
    ii = ii+1;
end 
x_monte_carlo_A_G_E = x(end);

x_monte_accuracy_g = x;
y_monte_accuracy_g = y;
u_monte_accuracy_g = u;
v_monte_accuracy_g = v;

figure(3)
subplot(2,2,3)
plot(x_monte_g,y_monte_g,'g')
hold on;
plot(x_monte_accuracy_g,y_monte_accuracy_g,'k-')
%plot(x_c,y_c,'*k')
hold off
title("Monte-Carlo accuracy [Grain]");
xlabel("X distance [m]");
ylabel("Y distance [m]");
%xlim([0.49,0.56])
legend("Grain Particles","Accuracy of Grain Particles");

clear x y u v
 
%  ----------------Monte-carlo accuracy [chaff]-----------------------%
rho_p=rho_c ;
d_p =d_c ;
N =100;
x=x0;
y=y0;
u=up0;
v=vp0;
y_c = -0.5;
dt=1e-4;
ii =1;

while(y(ii)>=y_c)
    [x(ii+1),y(ii+1),u(ii+1),v(ii+1)] = monte_carlo(x(ii),y(ii),u(ii),v(ii));
    ii = ii+1;
end 
x_monte_carlo_A_C_E = x(end);

x_monte_accuracy_c = x;
y_monte_accuracy_c = y;
u_monte_accuracy_c = u;
v_monte_accuracy_c = v;

figure(3)
subplot(2,2,4)
plot(x_monte_c,y_monte_c,'r')
hold on;
plot(x_monte_accuracy_c,y_monte_accuracy_c,'k-')
plot(x_c,y_c,'*k')
hold off
title('Monte-Carlo accuracy [chaff]');
xlabel('X-distance [m]')
ylabel('Y-diatance [m]')
xlim([0.49,0.64])
legend("Chaff Particles","Accuracy of Chaff Particles",'separation point');
% To increase the accuracy of the method, we use this 
set(gcf,'position',[10 10 1500 1500])
% to save the file we use, 
saveas(gcf,'Accuracy[G&C] with Euler & MC .jpg')
clear x y u v

%-------------------------Order of accuracy---------------------------%
%[interpolation is used to fill-in missing data, smooth existing data, make predictions.]
disp('Task-2 Quatifying error')
%interpolation of euler method to find the exact xp @ y_c
euler_grain = interp1(y_euler_g(end-1:end),x_euler_g(end-1:end),y_c);
euleraccuracy_grain = interp1(y_euler_accuracy_g(end-1:end),x_euler_accuracy_g(end-1:end),y_c);
euler_chaff = interp1(y_euler_c(end-1:end),x_euler_c(end-1:end),y_c);
euleraccuracy_chaff  = interp1(y_euler_accuracy_c(end-1:end),x_euler_accuracy_c(end-1:end),y_c);

%interpolation of Monte_carlo to find the exact xp @ y_c
monte_carlo_grain = interp1(y_monte_g(end-1:end),x_monte_g(end-1:end),y_c);
montecarlo_accuracy_grain = interp1(y_monte_accuracy_g(end-1:end),x_monte_accuracy_g(end-1:end),y_c);
monte_carlo_chaff = interp1(y_monte_c(end-1:end),x_monte_c(end-1:end),y_c);
montecarlo_accuracy_chaff = interp1(y_monte_accuracy_c(end-1:end),x_monte_accuracy_c(end-1:end),y_c);

% Error estimation
error_euler_grain = abs((euler_grain-euleraccuracy_grain)/euleraccuracy_grain);
error_euler_chaff = abs((euler_chaff-euleraccuracy_chaff)/euleraccuracy_chaff);
error_montecarlo_grain = abs((monte_carlo_grain-montecarlo_accuracy_grain)/montecarlo_accuracy_grain);
error_montecarlo_chaff = abs((monte_carlo_chaff-montecarlo_accuracy_chaff)/montecarlo_accuracy_chaff);

% error percentage
Absolute_error_euler_grain = error_euler_grain*100
Absolute_error_euler_chaff = error_euler_chaff*100
Absolute_error_montecarlo_grain = error_montecarlo_grain*100
Absolute_error_montecarlo_chaff = error_montecarlo_chaff*100

 %-------------------------Relative error estimation----------------------%
%--------------------------------Grains-----------------------------------%
rho_p=rho_g ;
d_p =d_g ;
N =[10,100,1000];
dt=1e-2;
for jj = 1:length(N)
   ii=1;
   x=x0;
   y=y0;
   u=up0;
   v=vp0;
   while(y(ii)>=y_c)
        [x(ii+1),y(ii+1),u(ii+1),v(ii+1)] =log_error(x(ii),y(ii),u(ii),v(ii),dt,N(jj));
        ii = ii+1;
   end 
   x_end_MC_error_G(jj)= x(end); 
   error_gr(jj) = abs(x_monte_carlo_A_G_E - x_end_MC_error_G(jj))/(x_monte_carlo_A_G_E);
end
figure(4)
subplot(1,2,1)
loglog(N,error_gr,'b--')
title('Study of the accuracy of the Monte-Carlo method[grain]@dt=0.001')
xlabel('Number of samples[N]')
ylabel('Relative error')
grid on

clear x y u v
%-------------------------------chaff------------------------------------%
rho_p=rho_c ;
d_p =d_c ;
N =[10,100,1000];
dt=1e-2;
for jj = 1:length(N)
   ii=1;
   x=x0;
   y=y0;
   u=up0;
   v=vp0;
   while(y(ii)>=y_c)
        [x(ii+1),y(ii+1),u(ii+1),v(ii+1)] =log_error(x(ii),y(ii),u(ii),v(ii),dt,N(jj));
        ii = ii+1;
   end 
   x_end_MC_error_C(jj)= x(end); 
   error_ch(jj) = abs(x_monte_carlo_A_C_E - x_end_MC_error_C(jj))/(x_monte_carlo_A_C_E);
end
figure(4)
subplot(1,2,2)
loglog(N,error_ch,'b--')
title('Study of the accuracy of the Monte-Carlo method[chaff]@dt=0.001')
xlabel('Number of samples[N]')
ylabel('Relative error')
grid on
set(gcf,'position',[10 10 1500 1500])
saveas(gcf,'loglog plot with MC.jpg')
clear x y u v

 %---------------------------Case[3]---------------------------------%
            %Euler approach is used to solve this task%

% Derive terminal velocity of each particle, @Rep>800, Newton regime(turbulent)(Cd=0.44) ) 

%Grains falling from angle_2 = -85, falls at the furtheset distance x1(end)for grains
%Chaffs falling from angle_1 = -95, falls at the nearest distance x2(end) for chaffs
% when(x1(end))<(x2(end)) =====> complete separation.

%---------------------------Grain-------------------------------------%
                            %[a]
times =10;
NOS = 10;  %number of grain sample 
grain_count_bin2 = 0;
for k =1:times
    for jj = 1:NOS
        rho_p = rho_g ;% the density of grain particle remains constant throghout the case[3] 
        sigma_d_g = 1 * 10^-3;
        % The normrand function is used to generate the normal distributon of grain diameter for given sigma_d_g.
        rand_d_g = normrnd(d_g,sigma_d_g);
        d_p = rand_d_g ;
        angle_gr = deg2rad((-95)+10*rand());
        % Intially, we assume that the particle fall in turbulent regime because the system equation for solving laminar case is non linear and it makes the calcualtion trivial 
        Ut_g_estimate  = sqrt(4*d_p*(rho_f-rho_p)*g/(3*rho_f*0.44)); %grain terminal velocity estimate.
        % with this Ut_estimate, we can find the actual regime in which the particle will flow 
        Re_p= rho_f* Ut_g_estimate * d_p/vis_f;
        if Re_p < 800
            Cd = (24/Re_p)*(1+(0.15*(Re_p)^(0.687)));
        else 
            Cd = 0.44;
        end
        Ut_g = sqrt(4*d_p*(rho_f-rho_p)*g/(3*rho_f*Cd));
        x=x0;
        y=y0;
        u= Ut_g*cos(angle_gr);
        v= Ut_g*sin(angle_gr);
        dt=1e-2;
        if (d_p > 0)
            ii=1;
            while (y(ii)>=y_c)
                [du1,dv1] = particle_velocity_ME(x(ii),y(ii),u(ii),v(ii));
                x(ii+1)=x(ii)+ dt*(u(ii));
                y(ii+1)=y(ii)+ dt*(v(ii));
                u(ii+1) = u(ii)+dt*(du1);
                v(ii+1)= v(ii)+dt*(dv1);
                ii= ii+1;
            end
            x = x(1:ii);
            y = y(1:ii);
            x_end_gr(jj)= x(end);
            if (x_end_gr(jj) >= x_c)
            grain_count_bin2 = grain_count_bin2 + 1;
            end 
            figure(5)
            plot(x,y,'r')
            hold on
        end 
    end
end 

grain_count_bin2 = grain_count_bin2/times;
proportion_grain = (grain_count_bin2 / NOS)*100;
fprintf('Proportions of Grain fallen into bin-2 is: %f%% \n\n',proportion_grain);


     %--------------------------------chaff-----------------------------%
                                      %[b]
times =10;
NOS = 10;  %number of chaff sample 
% The chaff particle is uniformly distributed in between the diameter range 2.0mm to 3.0mm
d_c_1 = 2*10^-3;
d_c_2 = 5*10^-3;
grain_count_bin1 = 0;
for k =1:times
    for jj = 1:NOS
        sigma_rho_c = 20;%[kg/m^3]
        rand_rho_c = normrnd(rho_c,sigma_rho_c);
        rho_p = rand_rho_c ;% the density of chaff particle is normally distributed.
        % The unifrnd function is used to generate the uniform distributon of grain diameter for given range
        rand_d_c = unifrnd(d_c_1,d_c_2);
        d_p = rand_d_c ;
        Ut_c_estimate = sqrt(4*d_p*(rho_f-rho_p)*g/(3*rho_f*0.44)); %chaff
        angle_ch = deg2rad((-95)+10*rand());
        Re_p= rho_f* Ut_c_estimate * d_p/vis_f;
        if Re_p < 800
            Cd = (24/Re_p)*(1+(0.15*(Re_p)^(0.687)));
        else 
            Cd = 0.44;
        end
        Ut_c = sqrt(4*d_p*(rho_f-rho_p)*g/(3*rho_f*Cd));
        x=x0;
        y=y0;
        u= Ut_g*cos(angle_ch);
        v= Ut_g*sin(angle_ch);
        dt=1e-2;        
        if (rho_p > rho_f) &&(d_p > 0)
            ii=1;
            while (y(ii)>=y_c)
                [du1,dv1] = particle_velocity_ME(x(ii),y(ii),u(ii),v(ii));
                x(ii+1)=x(ii)+ dt*(u(ii));
                y(ii+1)=y(ii)+ dt*(v(ii));
                u(ii+1) = u(ii)+dt*(du1);
                v(ii+1)= v(ii)+dt*(dv1);
                ii= ii+1;
            end
            x = x(1:ii);
            y = y(1:ii);
            x_end_ch(jj)= x(end);
            if (x_end_ch(jj) <= x_c)
            grain_count_bin1 = grain_count_bin1 + 1;
            end 
            figure(5)
            plot(x,y,'b')
            hold on
            xlabel('x distance [m]');
            ylabel('y distance [m]');
        end 
    end
end 
figure(5)
hold off
title('particle trajectory for varying diameter & density')
legend('grain particle','chaff particle')
set(gcf,'position',[10 10 1500 1500])

grain_count_bin1 = grain_count_bin1/times;
proportion_chaff = (grain_count_bin1 / NOS)*100;
fprintf('Proportions of chaff fallen into bin-1 is: %f%% \n\n',proportion_chaff);

%--------------------------------optimal distance of x_c---------------%
                                    %[c]
%clear all
x_minimum = 0.53;
x_maximum = 0.58;
step = 0.001;
x_c = x_minimum:step:x_maximum;
grain_count_bin2=zeros(1,length(x_c));  
grain_count_bin1=zeros(1,length(x_c));  
for i=1:length(x_c)
    NOS = 100;  %number of grain sample
    for jj = 1:NOS
        rho_p = rho_g ;% the density of grain particle remains constant throghout the case[3]
        sigma_d_g = 1 * 10^-3;
        % The normrand function is used to generate the normal distributon of grain diameter for given sigma_d_g.
        rand_d_g = normrnd(d_g,sigma_d_g);
        d_p = rand_d_g ;
        angle_gr = deg2rad((-95)+10*rand());
        % Intially, we assume that the particle fall in turbulent regime because the system equation for solving laminar case is non linear and it makes the calcualtion trivial 
        Ut_g_estimate  = sqrt(4*d_p*(rho_f-rho_p)*g/(3*rho_f*0.44)); %grain terminal velocity estimate.
        % with this Ut_estimate, we can find the actual regime in which the particle will flow 
        Re_p= rho_f* Ut_g_estimate * d_p/vis_f;
        if Re_p < 800
            Cd = (24/Re_p)*(1+(0.15*(Re_p)^(0.687)));
        else 
            Cd = 0.44;
        end
        Ut_g = sqrt(4*d_p*(rho_f-rho_p)*g/(3*rho_f*Cd));
        x=x0;
        y=y0;
        u=Ut_g*cos(angle_gr);
        v=Ut_g*sin(angle_gr);
        dt=1e-2;
        if (d_p > 0)
            ii=1;
            while (y(ii)>=y_c)
                [du1,dv1] = particle_velocity_ME(x(ii),y(ii),u(ii),v(ii));
                x(ii+1)=x(ii)+ dt*(u(ii));
                y(ii+1)=y(ii)+ dt*(v(ii));
                u(ii+1) = u(ii)+dt*(du1);
                v(ii+1)= v(ii)+dt*(dv1);
                ii= ii+1;
            end
            x = x(1:ii);
            y = y(1:ii);
            x_end_gr(jj)= x(end);
            if (x_end_gr(jj) >= x_c(i))
                grain_count_bin2(i) = grain_count_bin2(i) + 1;
            end
        end
    end
    %------------------------chaff-----------------------%
    NOS = 100;  %number of chaff sample
    % The chaff particle is uniformly distributed in between the diameter range 2.0mm to 3.0mm
    d_c_1 = 2*10^-3;
    d_c_2 = 5*10^-3;
    for jj = 1:NOS
        sigma_rho_c = 20;%[kg/m^3]
        rand_rho_c = normrnd(rho_c,sigma_rho_c);
        rho_p = rand_rho_c ;% the density of chaff particle is normally distributed.
        % The unifrnd function is used to generate the uniform distributon of grain diameter for given range
        rand_d_c = unifrnd(d_c_1,d_c_2);
        d_p = rand_d_c ;
        Ut_c_estimate = sqrt(4*d_p*(rho_f-rho_p)*g/(3*rho_f*0.44)); %chaff
        angle_ch = deg2rad((-95)+10*rand());
        Re_p= rho_f* Ut_c_estimate * d_p/vis_f;
        if Re_p < 800
            Cd = (24/Re_p)*(1+(0.15*(Re_p)^(0.687)));
        else 
            Cd = 0.44;
        end
        Ut_c = sqrt(4*d_p*(rho_f-rho_p)*g/(3*rho_f*Cd));
        x=x0;
        y=y0;
        u= Ut_c*cos(angle_ch);
        v= Ut_c*sin(angle_ch);
        dt=1e-2;
        if (rho_p > rho_f) &&(d_p > 0) % check the randomization gives a practical value 
            ii=1;
            while (y(ii)>=y_c)
                [du1,dv1] = particle_velocity_ME(x(ii),y(ii),u(ii),v(ii));
                x(ii+1)=x(ii)+ dt*(u(ii));
                y(ii+1)=y(ii)+ dt*(v(ii));
                u(ii+1) = u(ii)+dt*(du1);
                v(ii+1)= v(ii)+dt*(dv1);
                ii= ii+1;
            end
            x = x(1:ii);
            y = y(1:ii);
            x_end_ch(jj)= x(end);
            if (x_end_ch(jj) <= x_c(i))
                grain_count_bin1(i) = grain_count_bin1(i) + 1; 
            end
        end
    end   
end
% To check the proportion of particle in the irrespective bin 
    proportion_chaff= (grain_count_bin1 / NOS);
    proportion_grain = (grain_count_bin2 / NOS);
figure(6)
plot(x_c,proportion_grain,'b',x_c,proportion_chaff,'r')
title('Proportional distribution of both Grain and chaff ')
xlabel('x distance [m]')
ylabel('y distance [m]')
legend('Grain proportion','chaff proportion')
saveas(gcf,'probability plot.jpg')

toc
computation_run_time = toc; 

%--------------------------- functions --------------------------------%
 %-----------------function for fluid velcity ------------------------%

function [F_tu,F_tv]= particle_velocity_ME(x,y,up,vp)  
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
F_drag_y = (0.5*pi*rho_f*(d_p^2)*Cd*sqrt((uf-up)^2 + (vf-vp)^2)*(vf-vp));

F_tu = F_drag_x/(rho_p*V_p);
F_tv = ( F_gravity + F_buoyancy+F_drag_y )/(rho_p*V_p);
end

%----------------Function for monte carlo simulation---------------------%
function[x_end,y_end,u_end,v_end] = monte_carlo(x_ini,y_ini,u_ini,v_ini)
global dt N
dx = u_ini;
dy = v_ini;
t(1)=0;
t_rand = t(1)+ dt*rand(1,N);
% the t_rand will create a random set of values for both x and y 
x_rand =  x_ini+ dx*t_rand;
y_rand =  y_ini+ dy*t_rand;
for k =1:N
    [du(k),dv(k)] = particle_velocity_ME(x_rand(k),y_rand(k),u_ini,v_ini);
end 
mean_du = mean(du);
mean_dv = mean(dv);
x_end = x_ini + dx*dt;
y_end = y_ini + dy*dt;
u_end = u_ini + (mean_du)*dt;
v_end = v_ini + (mean_dv)*dt;
end 

    %-------------function for loglog ploting---------------------%
function[x_end,y_end,u_end,v_end] = log_error(x_ini,y_ini,u_ini,v_ini,dt,N)
dx = u_ini;
dy = v_ini;
t(1)=0;
t_rand = t(1)+ dt*rand(1,N);
% the t_rand will create a random set of values for both x and y 
x_rand =  x_ini+ dx*t_rand;
y_rand =  y_ini+ dy*t_rand;
for k =1:N
    [du(k),dv(k)] = particle_velocity_ME(x_rand(k),y_rand(k),u_ini,v_ini);
end 
mean_du = mean(du);
mean_dv = mean(dv);
x_end = x_ini + dx*dt;
y_end = y_ini + dy*dt;
u_end = u_ini + (mean_du)*dt;
v_end = v_ini + (mean_dv)*dt;
end 