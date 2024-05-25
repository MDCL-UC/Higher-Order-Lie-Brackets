clc
params

% Parameters for integral
P.k1 = k1;
P.k2 = k2;
P.k3 = k3; 
P.k4 = k4;


% Parameters for nu
P.p1 = p1; 
P.p2 = p2;  
P.p3 = p3; 
P.p4 = p4;
P.w  = w_val; 

% Parameters for vector fields
P.x_star = 1;  
P.h      = 5;  
P.alpha  = 1;  
P.H      = 1/5;  

% Parameters for simulation
P.simruntime = 40;  

% Initial values
I.x0 = 3;  
I.v0 = 0; 



%% Solver

% Gradient based ESC with 4 control inputs
[t_x_4i,X_4i] = ode45(@(t,X) gradientBasedESC4input(t,X,P), [0 P.simruntime], [I.x0;I.v0]);

% Gradient based ESC with 2 control inputs
[t_x_2i,X_2i] = ode45(@(t,X) gradientBasedESC2input(t,X,P), [0 P.simruntime], [I.x0;I.v0]);


% First order LBS (r = 2)
load nu_val2.mat
P.nu12 = nu_val2(12);
[t_xlbs1,X_lbs1] = ode45(@(t,X) gradientBasedFirstOrderLBS(X,P), [0 P.simruntime], [I.x0;I.v0]);

% Second order LBS (r = 3)
load nu_val3.mat
P.nu_val3 = double(values(nu_val3));
[t_xlbs2,X_lbs2] = ode45(@(t,X) gradientBasedSecondOrderLBS(X,P), [0 P.simruntime], [I.x0;I.v0]);
% save("gradientBasedSecondOrderLBS_4Controls_4thpowercostfunction","t_xlbs2","X_lbs2")

% Third order LBS (r = 4)
load beta_val.mat
P.all_beta = double(values(beta_val));
[t_xlbs3,X_lbs3] = ode45(@(t,X) gradientBasedThirdOrderLBS(X,P), [0 P.simruntime], [I.x0;I.v0]);

% Newton based ESC
load NewtonBased2019_4thordercost

%% Plots
fig = figure(1)
subplot(2,1,1)
plot(t_x_4i,X_4i(:,1),'g',"linewidth",2)
hold on
plot(t_xlbs1,X_lbs1(:,1),'r--',"linewidth",2)
plot(t_xlbs3,X_lbs3(:,1),'k--',"linewidth",2)
plot(t_Newton_based_ESC,x1_Newton_based_ESC,"m","linewidth",2)
plot(t_Newton_based_LBS,x1_Newton_based_LBS,"b--","linewidth",2)
title("State x vs time")
grid on
legend("Proposed ESC (m=4)","first-orderLBS (r = 2)","third-orderLBS(r = 4)","ESC in [1] (m=3)","LBS in [1] (r = 3)")
ylim([-1,8])
xlabel("Time")
ylabel("States")


subplot(2,1,2)
control_effort



%% Functions
function out = costFun(x,P)
x_star = P.x_star;
out = -P.H*(x-x_star)^4;  % <<<<<----- INPUT THIS
end



%% ESCs
function out = gradientBasedESC2input(t,X,P)

% Parameters
h = P.h;
w = P.w;
p1 = P.p1;
p2 = P.p2;
p3 = P.p3;
p4 = P.p4;
alpha = P.alpha;
k1 = P.k1;
k2 = P.k2;
k3 = P.k3;
k4 = P.k4;


% Variables
x = X(1);
v = X(2);

% Cost function
F = costFun(x,P);

% Dynamics
x_dot = (w^p1*(F-v)*sin(k1*w*t)+w^p2*alpha*cos(k2*w*t));
v_dot = h*(F-v);

% Output
out = [x_dot;v_dot];
end

function out = gradientBasedESC3input(t,X,P)

% Parameters
h = P.h;
w = P.w;
p1 = P.p1;
p2 = P.p2;
p3 = P.p3;
p4 = P.p4;
alpha = P.alpha;
k1 = P.k1;
k2 = P.k2;
k3 = P.k3;
k4 = P.k4;


% Variables
x = X(1);
v = X(2);

% Cost function
F = costFun(x,P);

% Dynamics
x_dot = (w^p1*(F-v)*sin(k1*w*t)+w^p2*alpha*cos(k2*w*t))+w^p3*alpha*sin(k3*w*t);
v_dot = h*(F-v);

% Output
out = [x_dot;v_dot];
end

function out = gradientBasedESC4input(t,X,P)

% Parameters
h = P.h;
w = P.w;
p1 = P.p1;
p2 = P.p2;
p3 = P.p3;
p4 = P.p4;
alpha = P.alpha;
k1 = P.k1;
k2 = P.k2;
k3 = P.k3;
k4 = P.k4;


% Variables
x = X(1);
v = X(2);

% Cost function
F = costFun(x,P);

% Dynamics
x_dot = (w^p1*(F-v)*sin(k1*w*t)+w^p2*alpha*cos(k2*w*t))+ w^p3*alpha*sin(k3*w*t)+ w^p4*alpha*cos(k4*w*t);
v_dot = h*(F-v);

% Output
out = [x_dot;v_dot];
end

%% LBSs
function out = gradientBasedFirstOrderLBS(X,P)
% Parameters
h = P.h;
a = P.alpha;
x_star = P.x_star;
H = P.H;
% nu12 = P.nu12;
nu12 = P.nu12;

% Variables
x = X(1);
v = X(2);

% Cost function
F = costFun(x,P);

% Lie brackets
f1f2 = 4*H*a*(x - x_star)^3;  % <<<<<----- INPUT THIS

% Dynamics
x_dot = nu12*f1f2;
v_dot = h*(F-v);

% Output
out = [x_dot;v_dot];
end

function out = gradientBasedSecondOrderLBS(X,P)
% Parameters
h = P.h;
x_star = P.x_star;
a = P.alpha;
H = P.H;
x = X(1);
v = X(2);
F = costFun(x,P);

all_nu = P.nu_val3;

% Lie brackets
f1f2_f1 = 12*H*a*(v + H*(x - x_star)^4)*(x - x_star)^2 - 16*H^2*a*(x - x_star)^6;
f1f2_f2 =  -12*H*a^2*(x - x_star)^2;

f = dictionary();
f(121) = f1f2_f1;
f(122) = f1f2_f2;
f(123) = f1f2_f2;
f(124) = f1f2_f2;

f(131) = f1f2_f1;
f(132) = f1f2_f2;
f(133) = f1f2_f2;
f(134) = f1f2_f2;

f(141) = f1f2_f1;
f(142) = f1f2_f2;
f(143) = f1f2_f2;
f(144) = f1f2_f2;

% The negative sign is present since we need f1_f1f2 instead of f1f2_f1.
all_f = -f.values();


% Dynamics
xdot1 = gradientBasedFirstOrderLBS(X,P);
xdot1 = xdot1(1);
xdot2 = sum(all_nu.*all_f); 
xdot = xdot1+xdot2;

vdot = h*(F-v);

% Output
out = [xdot;vdot];
end

function out = gradientBasedThirdOrderLBS(X,P)
% Or first order LBS
h = P.h;
H =P.H;
alpha = P.alpha;
x_star = P.x_star;
a = alpha;
x = X(1);
v = X(2);
F = costFun(x,P);

f1f2f1f1 = 4*H*(16*H^2*a*(x - x_star)^6 - 12*H*a*(v + H*(x - x_star)^4)*(x - x_star)^2)*(x - x_star)^3 - (48*H^2*a*(x - x_star)^5 - 12*H*a*(2*x - 2*x_star)*(v + H*(x - x_star)^4))*(v + H*(x - x_star)^4);
f1f2f1f2 = a*(48*H^2*a*(x - x_star)^5 - 12*H*a*(2*x - 2*x_star)*(v + H*(x - x_star)^4));
f1f2f2f1 = 48*H^2*a^2*(x - x_star)^5 - 12*H*a^2*(2*x - 2*x_star)*(v + H*(x - x_star)^4);
f1f2f2f2 = 12*H*a^3*(2*x - 2*x_star);


% f1f2f1f1 = f1f2f1f1;
% f1f2f1f2 = f1f2f1f2;
f1f2f1f3 = f1f2f1f2;
f1f2f1f4 = f1f2f1f2;

% f1f2f2f1 = f1f2f1f1;
% f1f2f2f2 = f1f2f1f2;
f1f2f2f3 = f1f2f2f2;
f1f2f2f4 = f1f2f2f2;

f1f2f3f1 = f1f2f2f1;
f1f2f3f2 = f1f2f2f2;
f1f2f3f3 = f1f2f2f2;
f1f2f3f4 = f1f2f2f2;

f1f2f4f1 = f1f2f2f1;
f1f2f4f2 = f1f2f2f2;
f1f2f4f3 = f1f2f2f2;
f1f2f4f4 = f1f2f2f2;


% 13
f1f3f1f1 = f1f2f1f1;
f1f3f1f2 = f1f2f1f2;
f1f3f1f3 = f1f2f1f2;
f1f3f1f4 = f1f2f1f2;

f1f3f2f1 = f1f2f2f1;
f1f3f2f2 = f1f2f2f2;
f1f3f2f3 = f1f2f2f2;
f1f3f2f4 = f1f2f2f2;

f1f3f3f1 = f1f2f2f1;
f1f3f3f2 = f1f2f2f2;
f1f3f3f3 = f1f2f2f2;
f1f3f3f4 = f1f2f2f2;

f1f3f4f1 = f1f2f2f1;
f1f3f4f2 = f1f2f2f2;
f1f3f4f3 = f1f2f2f2;
f1f3f4f4 = f1f2f2f2;

%14
f1f4f1f1 = f1f2f1f1;
f1f4f1f2 = f1f2f1f2;
f1f4f1f3 = f1f2f1f2;
f1f4f1f4 = f1f2f1f2;

f1f4f2f1 = f1f2f2f1;
f1f4f2f2 = f1f2f2f2;
f1f4f2f3 = f1f2f2f2;
f1f4f2f4 = f1f2f2f2;

f1f4f3f1 = f1f2f2f1;
f1f4f3f2 = f1f2f2f2;
f1f4f3f3 = f1f2f2f2;
f1f4f3f4 = f1f2f2f2;

f1f4f4f1 = f1f2f2f1;
f1f4f4f2 = f1f2f2f2;
f1f4f4f3 = f1f2f2f2;
f1f4f4f4 = f1f2f2f2;

all_vec_field = [f1f2f1f1 ;f1f2f1f2 ;f1f2f1f3 ;f1f2f1f4 ;
f1f2f2f1 ;f1f2f2f2;f1f2f2f3;f1f2f2f4 ;
f1f2f3f1 ;f1f2f3f2 ;f1f2f3f3 ;f1f2f3f4 ;
f1f2f4f1 ;f1f2f4f2 ;f1f2f4f3 ;f1f2f4f4 ;
%
f1f3f1f1 ;f1f3f1f2 ;f1f3f1f3 ;f1f3f1f4 ;
f1f3f2f1 ;f1f3f2f2 ;f1f3f2f3 ;f1f3f2f4 ;
f1f3f3f1 ;f1f3f3f2 ;f1f3f3f3 ;f1f3f3f4 ;
f1f3f4f1 ;f1f3f4f2 ;f1f3f4f3 ;f1f3f4f4 ;
%
f1f4f1f1 ;f1f4f1f2 ;f1f4f1f3 ;f1f4f1f4 ;
f1f4f2f1 ;f1f4f2f2 ;f1f4f2f3 ;f1f4f2f4 ;
f1f4f3f1 ;f1f4f3f2 ;f1f4f3f3 ;f1f4f3f4 ;
f1f4f4f1 ;f1f4f4f2 ;f1f4f4f3 ;f1f4f4f4 
];

% xdot1 = dynamicsEg4_LBS(x,P);
xdot1and2 = gradientBasedSecondOrderLBS(X,P);
xdot1and2 = xdot1and2(1);
xdot3 =  sum(P.all_beta.*all_vec_field);

x_dot = xdot1and2+xdot3;
v_dot = h*(F-v);
out = [x_dot;v_dot];
end
