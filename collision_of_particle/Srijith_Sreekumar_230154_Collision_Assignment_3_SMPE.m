%--------------------------Assignment 3---------------------------------%
%------------------Srijith Sreekumar "230154" ----------------------------%

close all;
clc;
tic;
%%%%%%%% All the function are mentioned at end of the script file%%%%%%%%

  %-----------------All values are in SI units------------------------%
%----------------------------Given data---------------------------------%
% Initial position of the particles.
x1i=[0 0 0];          x2i=[1.1 1.1 0];
% Initial velocities of the particles.
v1i=[0 0 0];          v2i=[-1.0 -1.0 0];

% Diameter of the particles
d1=1.0;             d2=1.0;
% Masses of the particles
m1=0.04;            m2=0.04; %Mass in Kg

% Gravitational constant 
g = -9.81;                   % [m/s2]

%Relative position of the particles
x1_2 = x1i-x2i;
%Relative velocities of the particles
v1_2 = v1i-v2i;

%--------------------------case[1]---------------------------------%
disp("------------------------Hard-Sphere Algorithm---------------------")
disp("Case[1]:")

% The collision time of the particle for a Hard-sphere algorithm:
minimum_collision_time = collisiontime(x1i,x2i,v1i,v2i,d1,d2);
fprintf('Minimum collision time = %f s \n\n',minimum_collision_time)


%--------------------------case[2]---------------------------------%
disp("Case[2]:")
%--------------------------Task[a]---------------------------------% 
disp("Task [2a]")   % The point of contact before collision.
x1 = x1i;
x2 = x2i;
v1 = v1i;
v2 = v2i;
t = 0;
dt = 1e-4;
d = d1;
count = 0;
i = 0;
while count==0
    i=i+1;
    x1 = x1 + (v1*dt);
    x2 = x2 + (v2*dt);
    t=t+dt;
    x12=x1-x2;
    if norm(x12)<=d
        count=count+100;
    end
end
x1_end = sprintf('%3f ',x1);
x2_end = sprintf('%3f ',x2);
fprintf('The final position of particle-1 before collision is: [%s] m \n',x1_end)
fprintf('The final position of particle-2 before collision is: [%s] m \n\n',x2_end)
% Since we have only one moving body (particle-2). The unit vector (n)of the
% relative position will yeild the direction of exact contact location.

x1_2_N =  x1_2; % to get the positive value of the unit vector
%unit vector
n = x1_2_N/norm(x1_2_N);
% multiplying with the radius of the particle 1 , will provide the
% magnitude to this direction vector(unit vector)

len_n = (d2/2)*n;

% Adding the centroid of the particle(stationary) to the len_n will proide the contact point co-ordinate 

cont_pt = x2i + len_n ;
contact_point = sprintf('%3f ',cont_pt);
fprintf('The contact point Co-ordinate before collision : [%s]m \n\n',contact_point)

clear x1 x2 v1 v2 x1_end x2_end count contact_point x1_2_N 

%--------------------------Task[b]---------------------------------% 
disp("Task [2b]") %the post-collision velocities of the two particles

%e - coefficient of restitution
e = 1.0;
[v1_final_c1,v2_fianl_c1] = finalvelocity_aftercollision(m1,m2,e,x1i,x2i,v1i,v2i);
disp ("For the coefficient of restitution: e = 1.0 ")
V1_c1 = sprintf('%3f ',v1_final_c1);
V2_c1 = sprintf('%3f ',v2_fianl_c1);
fprintf('The final velocity of particle-1 is: [%s] m/s \n',V1_c1)
fprintf('The final velocity of particle-2 is: [%s] m/s \n\n',V2_c1)

clear e
%--------------------------Task[c]---------------------------------% 
disp("Task [2c]")
e =0.75;
[v1_final_c2,v2_fianl_c2] = finalvelocity_aftercollision(m1,m2,e,x1i,x2i,v1i,v2i);
disp ("For the coefficient of restitution : e = 0.75 ")
V1_c2 = sprintf('%3f ',v1_final_c2);
V2_c2 = sprintf('%3f ',v2_fianl_c2);
fprintf('The final velocity of particle-1 is: [%s] m/s \n',V1_c2)
fprintf('The final velocity of particle-2 is: [%s] m/s \n\n',V2_c2)

clear e
%--------------------------Task[d]---------------------------------% 
disp("Task [2d]")
% Since in hard-sphere algorithm the collision is instaneous and binary.
% Mostly, the solution of task 2d is alomost 0.
dt=1e-5;
contact_time=0;
count =  0;
e=0.75;  
% just for reference collision time is drafted 
collision = collisiontime(x1i,x2i,v1i,v2i,d1,d2);
for t=0:dt:4
    if count==1  % for each time step, there is a possibility,thst it will run only once.
        mstar=m1*m2/(m1+m2); %reduced mass
        x1_2= x1i-x2i; %relative position
        n = x1_2/norm(x1_2); % collision normal
        v1_2=v1i-v2i; %relative inital velocity
        J =-(1+e)*dot(v1_2,n)*mstar.*n;
        v1f= v1i + (J/m);
        v2f= v2i - (J/m);
    end
    x1 = x1i + (v1i*dt);
    x2 = x2i + (v2i*dt);
    x1_2=x1-x2;
    if abs(norm(x1_2)-d1)<=0.01
        count=count+1;
        contact_time=contact_time+dt;
    end
end
fprintf('The contact time is: %.3f s\n\n',contact_time);

clear dt

%--------------------------Task[e]---------------------------------% 
disp("Task [2e]")
x1 = x1i;
x2 = x2i;
v1 = v1i;
v2 = v2i;
t = 0;
dt = 1e-4;
d = d1;
count = 0;
i = 0;
while count==0
    i=i+1;
    x1 = x1 + (v1*dt);
    x2 = x2 + (v2*dt);
    t=t+dt;
    x12=x1-x2;
    if norm(x12)<=d
        count=count+100;
    end
end
x1_end = sprintf('%3f ',x1);
x2_end = sprintf('%3f ',x2);
disp("Since in the hard-sphere algorithm the time of collision is binary and instantaneous.")
disp("The point after the collision will also be the same as the time before collision.")
fprintf('The final position of particle-1 after collision is: [%s] m \n',x1_end)
fprintf('The final position of particle-2 after collision is: [%s] m \n\n',x2_end)
% contact point coordinate
x1_2_N = - x1_2;
n = x1_2_N/norm(x1_2_N);
len_n = (d1/2)*n;
cont_pt = x1i + len_n ;
contact_point = sprintf('%3f ',cont_pt);
fprintf('The contact point Co-ordinate after collision : [%s]m \n\n\n',contact_point)

clear x1 x2 v1 v2 x1_end x2_end count contact_point x1_2_N 

%----------------------------case [3]-------------------------------%
fprintf("---------------------Soft-Sphere Algorithm---------------------\n\n")
fprintf("Case[3] - The function is mentioned at the bottom of the script file \n\n")
%----------------------------case [4]-------------------------------%
fprintf("Case[4] - The function is mentioned at the bottom of the script file \n\n")         
                         %Leapfrog algorithm

% The Leapfrog algorithm is applicable only from the point where the
% collision starts. That is the point when the colliding particle will excert
% some force on the stationary particle. This inturn would change the
% acceleation ,velocity and obviously the position of the particle would change.

%----------------------------case [5]-------------------------------%
%                       Soft-sphere algorithm                       %
disp("Case[5]")
%----------------------------Task [5a]-------------------------------%
disp("Task[5a]")
x1 = x1i;
x2 = x2i;
v1 = v1i;
v2 = v2i;
t = 0;
dt = 1e-4;
d = d1;
count = 0;
i = 0;
while count==0
    i=i+1;
    x1 = x1 + (v1*dt);
    x2 = x2 + (v2*dt);
    t=t+dt;
    x12=x1-x2;
    if norm(x12)<=d
        count=count+100;
    end
end
x1_end = sprintf('%3f ',x1);
x2_end = sprintf('%3f ',x2);
fprintf('The final position of particle-1 before collision is: [%s] m \n',x1_end)
fprintf('The final position of particle-2 before collision is: [%s] m \n\n',x2_end)
% contact point coordinate
x1_2_N = - x1_2;
n = x1_2_N/norm(x1_2_N);
len_n = (d1/2)*n;
cont_pt = x1i + len_n ;
contact_point = sprintf('%3f ',cont_pt);
fprintf('The contact point Co-ordinate before collision : [%s]m \n\n',contact_point)

clear x1 x2 v1 v2 x1_end x2_end count contact_point

%----------------------------Task [5b]-------------------------------%
disp("Task[5b]")
% For coefficient of resitution =1.0 , K_loading and k_unloading is same.
k = 750; 
e = 1.0;
cont_time = 0;
dt =  1e-4;
t_start = 0.3; % Since from case[1], we know the intial collision time of the particle will exactly start from 0.393 sec.
tmax = 2;% maximum duration
% Feeding the intial condition.
x1 = x1i;
x2 = x2i;
v1 = v1i;
v2 = v2i;
m = m1;
for ii = t_start:dt:tmax
    [x1,x2,v1,v2,~,~] = particle_motion(x1,x2,v1,v2,m,dt,k,d1,d2,e,cont_time);
end 
v1_end_c1 = sprintf('%3f ',v1);
v2_end_c1 = sprintf('%3f ',v2);
disp ("For the coefficient of restitution: e = 1.0 ")
fprintf('The final velocity of particle-1 is: [%s] m/s \n',v1_end_c1)
fprintf('The final velocity of particle-2 is: [%s] m/s \n\n',v2_end_c1)

clear k x1 x2 v1 v2 cont_time

%----------------------------Task [5c]-------------------------------%
disp("Task[5c]")

x1 = x1i;
x2 = x2i;
v1 = v1i;
v2 = v2i;
tmax = 2 ; % maximum duration
dt =  1e-4;
t_start = 0.3; % Since from case[1], we know the intial collision time of the particle will exactly start from 0.393 sec.
k = 750;
e = 0.75;
cont_time = 0;
for ii = t_start:dt:tmax
    [x1,x2,v1,v2,cont_time,~] = particle_motion(x1,x2,v1,v2,m,dt,k,d1,d2,e,cont_time);
end 
v1_end_c2 = sprintf('%3f ',v1);
v2_end_c2 = sprintf('%3f ',v2);
disp ("For the coefficient of restitution: e = 0.75")
fprintf('The final velocity of particle-1 is: [%s] m/s \n',v1_end_c2)
fprintf('The final velocity of particle-2 is: [%s] m/s \n\n',v2_end_c2)

%----------------------------Task [5d]-------------------------------%
disp("Task[5d]")
fprintf('The contact time for soft-sphere for e=0.75 is %f s \n\n',cont_time);

clear k x1 x2 v1 v2 cont_time 

%----------------------------Task [5e]-------------------------------%
disp("Task[5e]")

x1 = x1i;
x2 = x2i;
v1 = v1i;
v2 = v2i;
tmax = 1.5 ; % maximum duration
dt =  1e-4;
t_start = 0.3; % Since from case[1], we know the intial collision time of the particle will exactly start from 0.393 sec.
k = 750;
e = 1.0;
cont_time = 0;
contact_duration_i = 0;
contact_duration_f = 0;
for ii = t_start:dt:tmax
    [x1,x2,v1,v2,cont_time,overlap] = particle_motion(x1,x2,v1,v2,m,dt,k,d1,d2,e,cont_time);
    if overlap > 0
        contact_duration_i = 1; % the overlap is on. it stays true for the remainig period till loop ends 
    end
    % once the particle are no longer in contact. That particular is instant  has to be recorded  
    if overlap < 0 && contact_duration_i 
        contact_duration_f = contact_duration_f + 1 ; % this will increament untill the end of loop
    end 
    if contact_duration_f == 1 
        x12 = x2 - x1; % x1 ref to x2
        n = x12/norm(x12);
        len_n = (d1/2)*n;
        cont_pt = x1 + len_n ;
        x1_end_Ac = sprintf('%3f ',x1);
        x2_end_Ac = sprintf('%3f ',x2);
    end 
end 

fprintf('The final position of particle-1 after collision is: [%s] m \n',x1_end_Ac)
fprintf('The final position of particle-2 after collision is: [%s] m \n\n',x2_end_Ac)
contact_point = sprintf('%3f ',cont_pt);
fprintf('The contact point Co-ordinate after collision : [%s]m \n\n',contact_point)

% CPU run time
cpu_runtime = toc;
fprintf('The total CPU runtime is : %f s \n\n',cpu_runtime)



  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%Functions%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Function for finding the final velocities after collision 
function [v1_f,v2_f]=finalvelocity_aftercollision(m1,m2,e,x1_i,x2_i,v1_i,v2_i)
     mstar = (m1*m2)/(m1+m2); %Reduced mass
     x12 = x1_i-x2_i; %Relative position between the two particles
     n = x12/norm(x12); %Collision normal
     v12 = v1_i-v2_i; %Relative velocity between the two particles initially
     J = -(1+e)*dot(v12,n)*mstar*n; %Momentum exchange between particles
     %Individual final velocities of the particles
     v1_f = v1_i+(J/m1);
     v2_f = v2_i-(J/m2);
end

% Function for finding the minimum collision time 

                % reference Prof. code.
function mincollision = collisiontime(x1_i,x2_i,v1_i,v2_i,d1,d2)
    mincollision = -1;
    x12 = x1_i - x2_i;
    v12 = v1_i - v2_i;  %Relative velocities 
    b = dot(x12,v12);
    a = dot(v12,v12);
    c = dot(x12,x12) - (0.5*(d1+d2))^2;
    if (b>0) 
        return;
    end
    D = b^2 - a*c;
    if (D<0)
	    return;
    end
    mincollision = (-b - sqrt(D))/a;
end 

% Function for overlap and normal estimation(Case[3])
% Since tow arbitary location are being used to the result x1 and x2 values
% are used as an input aruguments 
function [F,overlap,normal,cont_time] = overlap_of_particle(x1,x2,v1,v2,k,d1,d2,e,cont_time,dt)
x12=x1-x2;
normal = x12/norm(x12);
delta = ((d1+d2)/2)-norm(x12);
v12 = v1-v2;
    %%%%%excecution on of force%%%%
    if delta > 0 % there is a collision .
          if dot(v12,normal)>0
              k = k*e^2;
              F = -k*delta*normal;
              cont_time = cont_time+dt; % updates the contact time of the particle 
          else 
               F = -k*delta*normal;
               cont_time = cont_time+dt;
          end 
    else 
        F = [0 0 0];        
    end 
overlap = delta;
end

% Functions for updation of position of the particles
% Velocity of the particles in two parts(predictor-corrector) using Leapfrog algorithm
function [x1,x2,v1,v2,cont_time,overlap] = particle_motion(x1,x2,v1,v2,m,dt,k,d1,d2,e,cont_time)
[F,~,~,~] = overlap_of_particle(x1,x2,v1,v2,k,d1,d2,e,cont_time,dt);
% for particle 1
v1= v1-((dt*F)/(2*m));    % Predictor velocity
x1= x1+v1*dt;             % Particle positionx_1=x+v*dt; % Particle position-1
% for particle 2 
v2= v2+((dt*F)/(2*m));    % Predictor velocity
x2= x2+v2*dt;             % Particle positionx_1=x+v*dt; % Particle position-2
%For the estimation of new accelaration @ a(i+1).
[F_new,overlap,~,cont_time] = overlap_of_particle(x1,x2,v1,v2,k,d1,d2,e,cont_time,dt); 
v1 = v1-((dt*F_new)/(2*m)); % Corrector velocity-1
v2 = v2+((dt*F_new)/(2*m)); % Corrector velocity-2
end
