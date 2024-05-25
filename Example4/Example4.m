clc
params
P.a = 1;
P.H = 1/3;

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

% Simulation
x0=[2];
P.simruntime = 20;

options=odeset('RelTol',1e-10,'Stats','on');
% ESC
[t_esc2grush,x_esc2grush]= ode15s(@(t,x) dynamicsEg4_ESC2Grush(t,x,P), [0, P.simruntime], x0,options);

% ESC
[t_esc2,x_esc2]= ode15s(@(t,x) dynamicsEg4_ESC2(t,x,P), [0, P.simruntime], x0,options);

% ESC
[t_esc4,x_esc4]= ode15s(@(t,x) dynamicsEg4_ESC4(t,x,P), [0, P.simruntime], x0, options);

% First order LBS Grush
load nu_val2.mat
[t_lbsgrush,x_lbsgrush]= ode15s(@(t,x) dynamicsEg4_LBS_Grush(x,P), [0, P.simruntime], x0,options);

% First order LBS corresponding to 4 control inputs
load nu_val2.mat
P.nu12 = nu_val2(12);
P.nu34 = nu_val2(34);
[t_lbs,x_lbs]= ode15s(@(t,x) dynamicsEg4_LBS(x,P), [0, P.simruntime], x0,options);

% First order LBS corresponding to 2 control inputs
[t_lbs_2control,x_lbs_2control]= ode15s(@(t,x) dynamicsEg4_LBS_2control(x,P), [0, P.simruntime], x0,options);

% Second order LBS
load nu_val3.mat
P.nu_val3 = double(values(nu_val3));
[t_lbs_2nd,x_lbs_2nd]=ode15s(@(t,x) dynamicsEg4_LBS2ndOrder(x,P), [0, P.simruntime], x0, options);


% Third order LBS
load beta_val.mat
all_beta = double(values(beta_val));
[t_lbs_3rd,x_lbs_3rd]=ode23(@(t,x) dynamicsEg4_LBS3rdOrder(x,P,all_beta), [0, P.simruntime], x0);

%% Plots
% Plots
figure(1)
subplot(2,1,1)
plot(t_esc2grush,x_esc2grush,'m',"LineWidth",2)
hold on
plot(t_lbsgrush,x_lbsgrush,'b--',"LineWidth",2)
plot(t_esc4,x_esc4,'g',"LineWidth",2)
plot(t_lbs,x_lbs,"color",[0.3010 0.7450 0.9330],"LineWidth",2)
plot(t_lbs_2nd,x_lbs_2nd,'r--',"LineWidth",2)
plot(t_lbs_3rd,x_lbs_3rd,'k--',"LineWidth",2)
grid on
legend("ESC in [22] (m=2)","LBS in [22] (r =2)", "Proposed ESC (m=4)", "LBS of Proposed ESC (r=2)", "LBS of Proposed ESC (r=3)", "LBS of Proposed ESC (r=4)")
xlabel("Time")
ylabel("States")
title("State vs time")


subplot(2,1,2)
control_effort
fontsize(10,"points")




%%
function out=dynamicsEg4_ESC2Grush(t,x,P)

k1 = P.k1;
k2 = P.k2;

p1 = P.p1;
p2 = P.p2;


H = P.H;
w = P.w;

J = H*(x - 1)^4;

phi = (1-exp(-(J)))/(1+exp(J));
psi = exp((J))+2*log(exp(J)-1);
u = sqrt(phi)*w^0.5*sin(psi)*cos(k1*w*t)+...
    sqrt(phi)*w^0.5*cos(psi)*sin(k2*w*t);

xdot = u;

out = xdot;

end

function out=dynamicsEg4_ESC2(t,x,P)

k1 = P.k1;
k2 = P.k2;

p1 = P.p1;
p2 = P.p2;


H = P.H;
w = P.w;

J = H*(x - 1)^4;

phi = (1-exp(-(J)))/(1+exp(J));
psi = exp((J))+2*log(exp(J)-1);
u = sqrt(phi)*w^p1*sin(psi)*cos(k1*w*t)+...
    sqrt(phi)*w^p2*cos(psi)*sin(k2*w*t);

xdot = u;

out = xdot;

end

function out=dynamicsEg4_ESC4(t,x,P)

k1 = P.k1;
k2 = P.k2;
k3 = P.k3;
k4 = P.k4;

p1 = P.p1;
p2 = P.p2;
p3 = P.p3;
p4 = P.p4;

H = P.H;
w = P.w;

J = H*(x - 1)^4;

phi = (1-exp(-(J)))/(1+exp(J));
psi = exp((J))+2*log(exp(J)-1);
u = sqrt(phi)*w^p1*sin(psi)*cos(k1*w*t)+...
    sqrt(phi)*w^p2*cos(psi)*sin(k2*w*t)+...
    sqrt(phi)*w^p3*sin(psi)*cos(k3*w*t)+...
    sqrt(phi)*w^p4*cos(psi)*sin(k4*w*t);

xdot = u;

out = xdot;

end

function out=dynamicsEg4_LBS_Grush(x,P)

% Parameters
H = P.H;

% Nu
nu12 = 0.5;

% LBs
f1f2 = -4*H*(x - 1)^3;

% Dynamics
xdot = nu12*f1f2;

% Outputs
out = xdot;
end

function out=dynamicsEg4_LBS_2control(x,P)
% Parameters
H = P.H;

% Nu
nu12 = double(P.nu12);



% LBs
f1f2 = -4*H*(x - 1)^3;

% Dynamics
xdot = nu12*f1f2;

% Outputs
out = xdot;
end

function out=dynamicsEg4_LBS(x,P)
% Parameters
H = P.H;

% Nu
nu12 = double(P.nu12);
nu34 = double(P.nu34);



% LBs
f1f2 = -4*H*(x - 1)^3;
f3f4 = f1f2;

% Dynamics
xdot = nu12*f1f2+nu34*f3f4;

% Outputs
out = xdot;
end

function out=dynamicsEg4_LBS2ndOrder(x,P)

x = x(1);
H = P.H;


% all_nu = double(values(nu));
all_nu = P.nu_val3;

% Vector field
f1f2_f1 = (12*H*exp((-0.5*H*(x - 1)^4))*sin(2*log(exp(H*(x - 1)^4) - 1) + exp(H*(x - 1)^4))*sqrt(exp(H*(x - 1)^4) - 1)*(x - 1)^2)/(exp(H*(x - 1)^4) + 1)^0.5 - 4*H*(x - 1)^3*((4*H*exp((0.5*H*(x - 1)^4))*cos(2*log(exp(H*(x - 1)^4) - 1) + exp(H*(x - 1)^4))*sqrt(exp(H*(x - 1)^4) + 1)*(x - 1)^3)/(exp(H*(x - 1)^4) - 1)^0.5 - (2*H*exp((0.5*H*(x - 1)^4))*sin(2*log(exp(H*(x - 1)^4) - 1) + exp(H*(x - 1)^4))*sqrt(exp(H*(x - 1)^4) - 1)*(x - 1)^3)/(exp(H*(x - 1)^4) + 1)^1.5 + (2*H*exp((-(0.5*H*(x - 1)^4)))*sin(2*log(exp(H*(x - 1)^4) - 1) + exp(H*(x - 1)^4))*(x - 1)^3)/((exp(H*(x - 1)^4) - 1)^0.5*(exp(H*(x - 1)^4) + 1)^0.5));
f1f2_f2 = 4*H*(x - 1)^3*((4*H*exp((0.5*H*(x - 1)^4))*sin(2*log(exp(H*(x - 1)^4) - 1) + exp(H*(x - 1)^4))*sqrt(exp(H*(x - 1)^4) + 1)*(x - 1)^3)/(exp(H*(x - 1)^4) - 1)^0.5 + (2*H*exp((0.5*H*(x - 1)^4))*cos(2*log(exp(H*(x - 1)^4) - 1) + exp(H*(x - 1)^4))*sqrt(exp(H*(x - 1)^4) - 1)*(x - 1)^3)/(exp(H*(x - 1)^4) + 1)^1 - (2*H*exp((-(0.5*H*(x - 1)^4)))*cos(2*log(exp(H*(x - 1)^4) - 1) + exp(H*(x - 1)^4))*(x - 1)^3)/((exp(H*(x - 1)^4) - 1)^0.5*(exp(H*(x - 1)^4) + 1)^0.5)) + (12*H*exp((-(0.5*H*(x - 1)^4)))*cos(2*log(exp(H*(x - 1)^4) - 1) + exp(H*(x - 1)^4))*sqrt(exp(H*(x - 1)^4) - 1)*(x - 1)^2)/(exp(H*(x - 1)^4) + 1)^0.5;

f = dictionary();
f(121) = f1f2_f1;
f(122) = f1f2_f2;
f(123) = f1f2_f1;
f(124) = f1f2_f2;

f(141) = f1f2_f1;
f(142) = f1f2_f2;
f(143) = f1f2_f1;
f(144) = f1f2_f2;

f(231)= f1f2_f1;
f(232) = f1f2_f2;
f(233) = f1f2_f1;
f(234) = f1f2_f2;

f(341) = f1f2_f1;
f(342) = f1f2_f2;
f(343)= f1f2_f1;
f(344) = f1f2_f2;


% The negative sign is present since we need f1_f1f2 instead of f1f2_f1.
all_f = -f.values();


% Dynamics
xdot1 = dynamicsEg4_LBS(x,P);
xdot2 = sum(all_nu.*all_f); 
xdot = xdot1+xdot2;

% Output
out = xdot;
end

function out = dynamicsEg4_LBS3rdOrder(x,P,all_beta)
H = P.H;




f1f2f1f1 = - (4*H*(x - 1)^3*((4*H*exp((0.5*H*(x - 1)^4))*cos(2*log(exp(H*(x - 1)^4) - 1) + exp(H*(x - 1)^4))*sqrt(exp(H*(x - 1)^4) + 1)*(x - 1)^3)/(exp(H*(x - 1)^4) - 1)^0.5 - (2*H*exp((0.5*H*(x - 1)^4))*sin(2*log(exp(H*(x - 1)^4) - 1) + exp(H*(x - 1)^4))*sqrt(exp(H*(x - 1)^4) - 1)*(x - 1)^3)/(exp(H*(x - 1)^4) + 1)^1.5 + (2*H*exp((-(0.5*H*(x - 1)^4)))*sin(2*log(exp(H*(x - 1)^4) - 1) + exp(H*(x - 1)^4))*(x - 1)^3)/((exp(H*(x - 1)^4) - 1)^0.5*(exp(H*(x - 1)^4) + 1)^0.5)) - (12*H*exp((-(0.5*H*(x - 1)^4)))*sin(2*log(exp(H*(x - 1)^4) - 1) + exp(H*(x - 1)^4))*sqrt(exp(H*(x - 1)^4) - 1)*(x - 1)^2)/(exp(H*(x - 1)^4) + 1)^0.5)*((4*H*exp((0.5*H*(x - 1)^4))*cos(2*log(exp(H*(x - 1)^4) - 1) + exp(H*(x - 1)^4))*sqrt(exp(H*(x - 1)^4) + 1)*(x - 1)^3)/(exp(H*(x - 1)^4) - 1)^0.5 - (2*H*exp((0.5*H*(x - 1)^4))*sin(2*log(exp(H*(x - 1)^4) - 1) + exp(H*(x - 1)^4))*sqrt(exp(H*(x - 1)^4) - 1)*(x - 1)^3)/(exp(H*(x - 1)^4) + 1)^1.5 + (2*H*exp((-(0.5*H*(x - 1)^4)))*sin(2*log(exp(H*(x - 1)^4) - 1) + exp(H*(x - 1)^4))*(x - 1)^3)/((exp(H*(x - 1)^4) - 1)^0.5*(exp(H*(x - 1)^4) + 1)^0.5)) - (exp((-(0.5*H*(x - 1)^4)))*sin(2*log(exp(H*(x - 1)^4) - 1) + exp(H*(x - 1)^4))*sqrt(exp(H*(x - 1)^4) - 1)*(4*H*(x - 1)^3*((4*H^2*exp(1.5*H*(x - 1)^4)*sin(2*log(exp(H*(x - 1)^4) - 1) + exp(H*(x - 1)^4))*(x - 1)^6)/((exp(H*(x - 1)^4) - 1)^0.5*(exp(H*(x - 1)^4) + 1)^1.5) - (12*H*exp((0.5*H*(x - 1)^4))*cos(2*log(exp(H*(x - 1)^4) - 1) + exp(H*(x - 1)^4))*sqrt(exp(H*(x - 1)^4) + 1)*(x - 1)^2)/(exp(H*(x - 1)^4) - 1)^0.5 + (16*H^2*exp(1.5*H*(x - 1)^4)*sin(2*log(exp(H*(x - 1)^4) - 1) + exp(H*(x - 1)^4))*(exp(H*(x - 1)^4) + 1)^1.5*(x - 1)^6)/(exp(H*(x - 1)^4) - 1)^1.5 + (6*H*exp((0.5*H*(x - 1)^4))*sin(2*log(exp(H*(x - 1)^4) - 1) + exp(H*(x - 1)^4))*sqrt(exp(H*(x - 1)^4) - 1)*(x - 1)^2)/(exp(H*(x - 1)^4) + 1)^1.5 - (6*H*exp((-(0.5*H*(x - 1)^4)))*sin(2*log(exp(H*(x - 1)^4) - 1) + exp(H*(x - 1)^4))*(x - 1)^2)/((exp(H*(x - 1)^4) - 1)^0.5*(exp(H*(x - 1)^4) + 1)^0.5) - (8*H^2*exp((0.5*H*(x - 1)^4))*cos(2*log(exp(H*(x - 1)^4) - 1) + exp(H*(x - 1)^4))*sqrt(exp(H*(x - 1)^4) + 1)*(x - 1)^6)/(exp(H*(x - 1)^4) - 1)^0.5 - (8*H^2*exp((0.5*H*(x - 1)^4))*cos(2*log(exp(H*(x - 1)^4) - 1) + exp(H*(x - 1)^4))*sqrt(exp(H*(x - 1)^4) + 1)*(x - 1)^6)/(exp(H*(x - 1)^4) - 1)^1.5 + (4*H^2*exp((0.5*H*(x - 1)^4))*sin(2*log(exp(H*(x - 1)^4) - 1) + exp(H*(x - 1)^4))*sqrt(exp(H*(x - 1)^4) - 1)*(x - 1)^6)/(exp(H*(x - 1)^4) + 1)^1.5 - (12*H^2*exp(1.5*H*(x - 1)^4)*sin(2*log(exp(H*(x - 1)^4) - 1) + exp(H*(x - 1)^4))*sqrt(exp(H*(x - 1)^4) - 1)*(x - 1)^6)/(exp(H*(x - 1)^4) + 1)^2.5 + (8*H^2*exp(1.5*H*(x - 1)^4)*cos(2*log(exp(H*(x - 1)^4) - 1) + exp(H*(x - 1)^4))*sqrt(exp(H*(x - 1)^4) + 1)*(x - 1)^6)/(exp(H*(x - 1)^4) - 1)^1.5 + (4*H^2*exp((-(0.5*H*(x - 1)^4)))*sin(2*log(exp(H*(x - 1)^4) - 1) + exp(H*(x - 1)^4))*(x - 1)^6)/((exp(H*(x - 1)^4) - 1)^0.5*(exp(H*(x - 1)^4) + 1)^0.5) + (4*H^2*exp((0.5*H*(x - 1)^4))*sin(2*log(exp(H*(x - 1)^4) - 1) + exp(H*(x - 1)^4))*(x - 1)^6)/((exp(H*(x - 1)^4) - 1)^0.5*(exp(H*(x - 1)^4) + 1)^1.5) + (4*H^2*exp((0.5*H*(x - 1)^4))*sin(2*log(exp(H*(x - 1)^4) - 1) + exp(H*(x - 1)^4))*(x - 1)^6)/((exp(H*(x - 1)^4) - 1)^1.5*(exp(H*(x - 1)^4) + 1)^0.5)) - 12*H*(x - 1)^2*((4*H*exp((0.5*H*(x - 1)^4))*cos(2*log(exp(H*(x - 1)^4) - 1) + exp(H*(x - 1)^4))*sqrt(exp(H*(x - 1)^4) + 1)*(x - 1)^3)/(exp(H*(x - 1)^4) - 1)^0.5 - (2*H*exp((0.5*H*(x - 1)^4))*sin(2*log(exp(H*(x - 1)^4) - 1) + exp(H*(x - 1)^4))*sqrt(exp(H*(x - 1)^4) - 1)*(x - 1)^3)/(exp(H*(x - 1)^4) + 1)^1.5 + (2*H*exp((-(0.5*H*(x - 1)^4)))*sin(2*log(exp(H*(x - 1)^4) - 1) + exp(H*(x - 1)^4))*(x - 1)^3)/((exp(H*(x - 1)^4) - 1)^0.5*(exp(H*(x - 1)^4) + 1)^0.5)) + (12*H*exp((-(0.5*H*(x - 1)^4)))*sin(2*log(exp(H*(x - 1)^4) - 1) + exp(H*(x - 1)^4))*(2*x - 2)*sqrt(exp(H*(x - 1)^4) - 1))/(exp(H*(x - 1)^4) + 1)^0.5 + (48*H^2*exp((0.5*H*(x - 1)^4))*cos(2*log(exp(H*(x - 1)^4) - 1) + exp(H*(x - 1)^4))*sqrt(exp(H*(x - 1)^4) + 1)*(x - 1)^5)/(exp(H*(x - 1)^4) - 1)^0.5 - (24*H^2*exp((-(0.5*H*(x - 1)^4)))*sin(2*log(exp(H*(x - 1)^4) - 1) + exp(H*(x - 1)^4))*sqrt(exp(H*(x - 1)^4) - 1)*(x - 1)^5)/(exp(H*(x - 1)^4) + 1)^0.5 - (24*H^2*exp((0.5*H*(x - 1)^4))*sin(2*log(exp(H*(x - 1)^4) - 1) + exp(H*(x - 1)^4))*sqrt(exp(H*(x - 1)^4) - 1)*(x - 1)^5)/(exp(H*(x - 1)^4) + 1)^1.5 + (24*H^2*exp((0.5*H*(x - 1)^4))*sin(2*log(exp(H*(x - 1)^4) - 1) + exp(H*(x - 1)^4))*(x - 1)^5)/((exp(H*(x - 1)^4) - 1)^0.5*(exp(H*(x - 1)^4) + 1)^0.5)))/(exp(H*(x - 1)^4) + 1)^0.5;
f1f2f1f2 = (4*H*(x - 1)^3*((4*H*exp((0.5*H*(x - 1)^4))*cos(2*log(exp(H*(x - 1)^4) - 1) + exp(H*(x - 1)^4))*sqrt(exp(H*(x - 1)^4) + 1)*(x - 1)^3)/(exp(H*(x - 1)^4) - 1)^0.5 - (2*H*exp((0.5*H*(x - 1)^4))*sin(2*log(exp(H*(x - 1)^4) - 1) + exp(H*(x - 1)^4))*sqrt(exp(H*(x - 1)^4) - 1)*(x - 1)^3)/(exp(H*(x - 1)^4) + 1)^1.5 + (2*H*exp((-(0.5*H*(x - 1)^4)))*sin(2*log(exp(H*(x - 1)^4) - 1) + exp(H*(x - 1)^4))*(x - 1)^3)/((exp(H*(x - 1)^4) - 1)^0.5*(exp(H*(x - 1)^4) + 1)^0.5)) - (12*H*exp((-(0.5*H*(x - 1)^4)))*sin(2*log(exp(H*(x - 1)^4) - 1) + exp(H*(x - 1)^4))*sqrt(exp(H*(x - 1)^4) - 1)*(x - 1)^2)/(exp(H*(x - 1)^4) + 1)^0.5)*((4*H*exp((0.5*H*(x - 1)^4))*sin(2*log(exp(H*(x - 1)^4) - 1) + exp(H*(x - 1)^4))*sqrt(exp(H*(x - 1)^4) + 1)*(x - 1)^3)/(exp(H*(x - 1)^4) - 1)^0.5 + (2*H*exp((0.5*H*(x - 1)^4))*cos(2*log(exp(H*(x - 1)^4) - 1) + exp(H*(x - 1)^4))*sqrt(exp(H*(x - 1)^4) - 1)*(x - 1)^3)/(exp(H*(x - 1)^4) + 1)^1.5 - (2*H*exp((-(0.5*H*(x - 1)^4)))*cos(2*log(exp(H*(x - 1)^4) - 1) + exp(H*(x - 1)^4))*(x - 1)^3)/((exp(H*(x - 1)^4) - 1)^0.5*(exp(H*(x - 1)^4) + 1)^0.5)) - (exp((-(0.5*H*(x - 1)^4)))*cos(2*log(exp(H*(x - 1)^4) - 1) + exp(H*(x - 1)^4))*sqrt(exp(H*(x - 1)^4) - 1)*(4*H*(x - 1)^3*((4*H^2*exp(1.5*H*(x - 1)^4)*sin(2*log(exp(H*(x - 1)^4) - 1) + exp(H*(x - 1)^4))*(x - 1)^6)/((exp(H*(x - 1)^4) - 1)^0.5*(exp(H*(x - 1)^4) + 1)^1.5) - (12*H*exp((0.5*H*(x - 1)^4))*cos(2*log(exp(H*(x - 1)^4) - 1) + exp(H*(x - 1)^4))*sqrt(exp(H*(x - 1)^4) + 1)*(x - 1)^2)/(exp(H*(x - 1)^4) - 1)^0.5 + (16*H^2*exp(1.5*H*(x - 1)^4)*sin(2*log(exp(H*(x - 1)^4) - 1) + exp(H*(x - 1)^4))*(exp(H*(x - 1)^4) + 1)^1.5*(x - 1)^6)/(exp(H*(x - 1)^4) - 1)^1.5 + (6*H*exp((0.5*H*(x - 1)^4))*sin(2*log(exp(H*(x - 1)^4) - 1) + exp(H*(x - 1)^4))*sqrt(exp(H*(x - 1)^4) - 1)*(x - 1)^2)/(exp(H*(x - 1)^4) + 1)^1.5 - (6*H*exp((-(0.5*H*(x - 1)^4)))*sin(2*log(exp(H*(x - 1)^4) - 1) + exp(H*(x - 1)^4))*(x - 1)^2)/((exp(H*(x - 1)^4) - 1)^0.5*(exp(H*(x - 1)^4) + 1)^0.5) - (8*H^2*exp((0.5*H*(x - 1)^4))*cos(2*log(exp(H*(x - 1)^4) - 1) + exp(H*(x - 1)^4))*sqrt(exp(H*(x - 1)^4) + 1)*(x - 1)^6)/(exp(H*(x - 1)^4) - 1)^0.5 - (8*H^2*exp((0.5*H*(x - 1)^4))*cos(2*log(exp(H*(x - 1)^4) - 1) + exp(H*(x - 1)^4))*sqrt(exp(H*(x - 1)^4) + 1)*(x - 1)^6)/(exp(H*(x - 1)^4) - 1)^1.5 + (4*H^2*exp((0.5*H*(x - 1)^4))*sin(2*log(exp(H*(x - 1)^4) - 1) + exp(H*(x - 1)^4))*sqrt(exp(H*(x - 1)^4) - 1)*(x - 1)^6)/(exp(H*(x - 1)^4) + 1)^1.5 - (12*H^2*exp(1.5*H*(x - 1)^4)*sin(2*log(exp(H*(x - 1)^4) - 1) + exp(H*(x - 1)^4))*sqrt(exp(H*(x - 1)^4) - 1)*(x - 1)^6)/(exp(H*(x - 1)^4) + 1)^2.5 + (8*H^2*exp(1.5*H*(x - 1)^4)*cos(2*log(exp(H*(x - 1)^4) - 1) + exp(H*(x - 1)^4))*sqrt(exp(H*(x - 1)^4) + 1)*(x - 1)^6)/(exp(H*(x - 1)^4) - 1)^1.5 + (4*H^2*exp((-(0.5*H*(x - 1)^4)))*sin(2*log(exp(H*(x - 1)^4) - 1) + exp(H*(x - 1)^4))*(x - 1)^6)/((exp(H*(x - 1)^4) - 1)^0.5*(exp(H*(x - 1)^4) + 1)^0.5) + (4*H^2*exp((0.5*H*(x - 1)^4))*sin(2*log(exp(H*(x - 1)^4) - 1) + exp(H*(x - 1)^4))*(x - 1)^6)/((exp(H*(x - 1)^4) - 1)^0.5*(exp(H*(x - 1)^4) + 1)^1.5) + (4*H^2*exp((0.5*H*(x - 1)^4))*sin(2*log(exp(H*(x - 1)^4) - 1) + exp(H*(x - 1)^4))*(x - 1)^6)/((exp(H*(x - 1)^4) - 1)^1.5*(exp(H*(x - 1)^4) + 1)^0.5)) - 12*H*(x - 1)^2*((4*H*exp((0.5*H*(x - 1)^4))*cos(2*log(exp(H*(x - 1)^4) - 1) + exp(H*(x - 1)^4))*sqrt(exp(H*(x - 1)^4) + 1)*(x - 1)^3)/(exp(H*(x - 1)^4) - 1)^0.5 - (2*H*exp((0.5*H*(x - 1)^4))*sin(2*log(exp(H*(x - 1)^4) - 1) + exp(H*(x - 1)^4))*sqrt(exp(H*(x - 1)^4) - 1)*(x - 1)^3)/(exp(H*(x - 1)^4) + 1)^1.5 + (2*H*exp((-(0.5*H*(x - 1)^4)))*sin(2*log(exp(H*(x - 1)^4) - 1) + exp(H*(x - 1)^4))*(x - 1)^3)/((exp(H*(x - 1)^4) - 1)^0.5*(exp(H*(x - 1)^4) + 1)^0.5)) + (12*H*exp((-(0.5*H*(x - 1)^4)))*sin(2*log(exp(H*(x - 1)^4) - 1) + exp(H*(x - 1)^4))*(2*x - 2)*sqrt(exp(H*(x - 1)^4) - 1))/(exp(H*(x - 1)^4) + 1)^0.5 + (48*H^2*exp((0.5*H*(x - 1)^4))*cos(2*log(exp(H*(x - 1)^4) - 1) + exp(H*(x - 1)^4))*sqrt(exp(H*(x - 1)^4) + 1)*(x - 1)^5)/(exp(H*(x - 1)^4) - 1)^0.5 - (24*H^2*exp((-(0.5*H*(x - 1)^4)))*sin(2*log(exp(H*(x - 1)^4) - 1) + exp(H*(x - 1)^4))*sqrt(exp(H*(x - 1)^4) - 1)*(x - 1)^5)/(exp(H*(x - 1)^4) + 1)^0.5 - (24*H^2*exp((0.5*H*(x - 1)^4))*sin(2*log(exp(H*(x - 1)^4) - 1) + exp(H*(x - 1)^4))*sqrt(exp(H*(x - 1)^4) - 1)*(x - 1)^5)/(exp(H*(x - 1)^4) + 1)^1.5 + (24*H^2*exp((0.5*H*(x - 1)^4))*sin(2*log(exp(H*(x - 1)^4) - 1) + exp(H*(x - 1)^4))*(x - 1)^5)/((exp(H*(x - 1)^4) - 1)^0.5*(exp(H*(x - 1)^4) + 1)^0.5)))/(exp(H*(x - 1)^4) + 1)^0.5;
f1f2f2f1 = (4*H*(x - 1)^3*((4*H*exp((0.5*H*(x - 1)^4))*sin(2*log(exp(H*(x - 1)^4) - 1) + exp(H*(x - 1)^4))*sqrt(exp(H*(x - 1)^4) + 1)*(x - 1)^3)/(exp(H*(x - 1)^4) - 1)^0.5 + (2*H*exp((0.5*H*(x - 1)^4))*cos(2*log(exp(H*(x - 1)^4) - 1) + exp(H*(x - 1)^4))*sqrt(exp(H*(x - 1)^4) - 1)*(x - 1)^3)/(exp(H*(x - 1)^4) + 1)^1.5 - (2*H*exp((-(0.5*H*(x - 1)^4)))*cos(2*log(exp(H*(x - 1)^4) - 1) + exp(H*(x - 1)^4))*(x - 1)^3)/((exp(H*(x - 1)^4) - 1)^0.5*(exp(H*(x - 1)^4) + 1)^0.5)) + (12*H*exp((-(0.5*H*(x - 1)^4)))*cos(2*log(exp(H*(x - 1)^4) - 1) + exp(H*(x - 1)^4))*sqrt(exp(H*(x - 1)^4) - 1)*(x - 1)^2)/(exp(H*(x - 1)^4) + 1)^0.5)*((4*H*exp((0.5*H*(x - 1)^4))*cos(2*log(exp(H*(x - 1)^4) - 1) + exp(H*(x - 1)^4))*sqrt(exp(H*(x - 1)^4) + 1)*(x - 1)^3)/(exp(H*(x - 1)^4) - 1)^0.5 - (2*H*exp((0.5*H*(x - 1)^4))*sin(2*log(exp(H*(x - 1)^4) - 1) + exp(H*(x - 1)^4))*sqrt(exp(H*(x - 1)^4) - 1)*(x - 1)^3)/(exp(H*(x - 1)^4) + 1)^1.5 + (2*H*exp((-(0.5*H*(x - 1)^4)))*sin(2*log(exp(H*(x - 1)^4) - 1) + exp(H*(x - 1)^4))*(x - 1)^3)/((exp(H*(x - 1)^4) - 1)^0.5*(exp(H*(x - 1)^4) + 1)^0.5)) - (exp((-(0.5*H*(x - 1)^4)))*sin(2*log(exp(H*(x - 1)^4) - 1) + exp(H*(x - 1)^4))*sqrt(exp(H*(x - 1)^4) - 1)*(12*H*(x - 1)^2*((4*H*exp((0.5*H*(x - 1)^4))*sin(2*log(exp(H*(x - 1)^4) - 1) + exp(H*(x - 1)^4))*sqrt(exp(H*(x - 1)^4) + 1)*(x - 1)^3)/(exp(H*(x - 1)^4) - 1)^0.5 + (2*H*exp((0.5*H*(x - 1)^4))*cos(2*log(exp(H*(x - 1)^4) - 1) + exp(H*(x - 1)^4))*sqrt(exp(H*(x - 1)^4) - 1)*(x - 1)^3)/(exp(H*(x - 1)^4) + 1)^1.5 - (2*H*exp((-(0.5*H*(x - 1)^4)))*cos(2*log(exp(H*(x - 1)^4) - 1) + exp(H*(x - 1)^4))*(x - 1)^3)/((exp(H*(x - 1)^4) - 1)^0.5*(exp(H*(x - 1)^4) + 1)^0.5)) + 4*H*(x - 1)^3*((16*H^2*exp(1.5*H*(x - 1)^4)*cos(2*log(exp(H*(x - 1)^4) - 1) + exp(H*(x - 1)^4))*(exp(H*(x - 1)^4) + 1)^1.5*(x - 1)^6)/(exp(H*(x - 1)^4) - 1)^1.5 + (12*H*exp((0.5*H*(x - 1)^4))*sin(2*log(exp(H*(x - 1)^4) - 1) + exp(H*(x - 1)^4))*sqrt(exp(H*(x - 1)^4) + 1)*(x - 1)^2)/(exp(H*(x - 1)^4) - 1)^0.5 + (6*H*exp((0.5*H*(x - 1)^4))*cos(2*log(exp(H*(x - 1)^4) - 1) + exp(H*(x - 1)^4))*sqrt(exp(H*(x - 1)^4) - 1)*(x - 1)^2)/(exp(H*(x - 1)^4) + 1)^1.5 - (6*H*exp((-(0.5*H*(x - 1)^4)))*cos(2*log(exp(H*(x - 1)^4) - 1) + exp(H*(x - 1)^4))*(x - 1)^2)/((exp(H*(x - 1)^4) - 1)^0.5*(exp(H*(x - 1)^4) + 1)^0.5) + (8*H^2*exp((0.5*H*(x - 1)^4))*sin(2*log(exp(H*(x - 1)^4) - 1) + exp(H*(x - 1)^4))*sqrt(exp(H*(x - 1)^4) + 1)*(x - 1)^6)/(exp(H*(x - 1)^4) - 1)^0.5 + (4*H^2*exp((0.5*H*(x - 1)^4))*cos(2*log(exp(H*(x - 1)^4) - 1) + exp(H*(x - 1)^4))*sqrt(exp(H*(x - 1)^4) - 1)*(x - 1)^6)/(exp(H*(x - 1)^4) + 1)^1.5 + (8*H^2*exp((0.5*H*(x - 1)^4))*sin(2*log(exp(H*(x - 1)^4) - 1) + exp(H*(x - 1)^4))*sqrt(exp(H*(x - 1)^4) + 1)*(x - 1)^6)/(exp(H*(x - 1)^4) - 1)^1.5 - (12*H^2*exp(1.5*H*(x - 1)^4)*cos(2*log(exp(H*(x - 1)^4) - 1) + exp(H*(x - 1)^4))*sqrt(exp(H*(x - 1)^4) - 1)*(x - 1)^6)/(exp(H*(x - 1)^4) + 1)^2.5 + (4*H^2*exp((-(0.5*H*(x - 1)^4)))*cos(2*log(exp(H*(x - 1)^4) - 1) + exp(H*(x - 1)^4))*(x - 1)^6)/((exp(H*(x - 1)^4) - 1)^0.5*(exp(H*(x - 1)^4) + 1)^0.5) - (8*H^2*exp(1.5*H*(x - 1)^4)*sin(2*log(exp(H*(x - 1)^4) - 1) + exp(H*(x - 1)^4))*sqrt(exp(H*(x - 1)^4) + 1)*(x - 1)^6)/(exp(H*(x - 1)^4) - 1)^1.5 + (4*H^2*exp((0.5*H*(x - 1)^4))*cos(2*log(exp(H*(x - 1)^4) - 1) + exp(H*(x - 1)^4))*(x - 1)^6)/((exp(H*(x - 1)^4) - 1)^0.5*(exp(H*(x - 1)^4) + 1)^1.5) + (4*H^2*exp((0.5*H*(x - 1)^4))*cos(2*log(exp(H*(x - 1)^4) - 1) + exp(H*(x - 1)^4))*(x - 1)^6)/((exp(H*(x - 1)^4) - 1)^1.5*(exp(H*(x - 1)^4) + 1)^0.5) + (4*H^2*exp(1.5*H*(x - 1)^4)*cos(2*log(exp(H*(x - 1)^4) - 1) + exp(H*(x - 1)^4))*(x - 1)^6)/((exp(H*(x - 1)^4) - 1)^0.5*(exp(H*(x - 1)^4) + 1)^1.5)) + (12*H*exp((-(0.5*H*(x - 1)^4)))*cos(2*log(exp(H*(x - 1)^4) - 1) + exp(H*(x - 1)^4))*(2*x - 2)*sqrt(exp(H*(x - 1)^4) - 1))/(exp(H*(x - 1)^4) + 1)^0.5 - (24*H^2*exp((-(0.5*H*(x - 1)^4)))*cos(2*log(exp(H*(x - 1)^4) - 1) + exp(H*(x - 1)^4))*sqrt(exp(H*(x - 1)^4) - 1)*(x - 1)^5)/(exp(H*(x - 1)^4) + 1)^0.5 - (48*H^2*exp((0.5*H*(x - 1)^4))*sin(2*log(exp(H*(x - 1)^4) - 1) + exp(H*(x - 1)^4))*sqrt(exp(H*(x - 1)^4) + 1)*(x - 1)^5)/(exp(H*(x - 1)^4) - 1)^0.5 - (24*H^2*exp((0.5*H*(x - 1)^4))*cos(2*log(exp(H*(x - 1)^4) - 1) + exp(H*(x - 1)^4))*sqrt(exp(H*(x - 1)^4) - 1)*(x - 1)^5)/(exp(H*(x - 1)^4) + 1)^1.5 + (24*H^2*exp((0.5*H*(x - 1)^4))*cos(2*log(exp(H*(x - 1)^4) - 1) + exp(H*(x - 1)^4))*(x - 1)^5)/((exp(H*(x - 1)^4) - 1)^0.5*(exp(H*(x - 1)^4) + 1)^0.5)))/(exp(H*(x - 1)^4) + 1)^0.5;
f1f2f2f2 = - (4*H*(x - 1)^3*((4*H*exp((0.5*H*(x - 1)^4))*sin(2*log(exp(H*(x - 1)^4) - 1) + exp(H*(x - 1)^4))*sqrt(exp(H*(x - 1)^4) + 1)*(x - 1)^3)/(exp(H*(x - 1)^4) - 1)^0.5 + (2*H*exp((0.5*H*(x - 1)^4))*cos(2*log(exp(H*(x - 1)^4) - 1) + exp(H*(x - 1)^4))*sqrt(exp(H*(x - 1)^4) - 1)*(x - 1)^3)/(exp(H*(x - 1)^4) + 1)^1.5 - (2*H*exp((-(0.5*H*(x - 1)^4)))*cos(2*log(exp(H*(x - 1)^4) - 1) + exp(H*(x - 1)^4))*(x - 1)^3)/((exp(H*(x - 1)^4) - 1)^0.5*(exp(H*(x - 1)^4) + 1)^0.5)) + (12*H*exp((-(0.5*H*(x - 1)^4)))*cos(2*log(exp(H*(x - 1)^4) - 1) + exp(H*(x - 1)^4))*sqrt(exp(H*(x - 1)^4) - 1)*(x - 1)^2)/(exp(H*(x - 1)^4) + 1)^0.5)*((4*H*exp((0.5*H*(x - 1)^4))*sin(2*log(exp(H*(x - 1)^4) - 1) + exp(H*(x - 1)^4))*sqrt(exp(H*(x - 1)^4) + 1)*(x - 1)^3)/(exp(H*(x - 1)^4) - 1)^0.5 + (2*H*exp((0.5*H*(x - 1)^4))*cos(2*log(exp(H*(x - 1)^4) - 1) + exp(H*(x - 1)^4))*sqrt(exp(H*(x - 1)^4) - 1)*(x - 1)^3)/(exp(H*(x - 1)^4) + 1)^1.5 - (2*H*exp((-(0.5*H*(x - 1)^4)))*cos(2*log(exp(H*(x - 1)^4) - 1) + exp(H*(x - 1)^4))*(x - 1)^3)/((exp(H*(x - 1)^4) - 1)^0.5*(exp(H*(x - 1)^4) + 1)^0.5)) - (exp((-(0.5*H*(x - 1)^4)))*cos(2*log(exp(H*(x - 1)^4) - 1) + exp(H*(x - 1)^4))*sqrt(exp(H*(x - 1)^4) - 1)*(12*H*(x - 1)^2*((4*H*exp((0.5*H*(x - 1)^4))*sin(2*log(exp(H*(x - 1)^4) - 1) + exp(H*(x - 1)^4))*sqrt(exp(H*(x - 1)^4) + 1)*(x - 1)^3)/(exp(H*(x - 1)^4) - 1)^0.5 + (2*H*exp((0.5*H*(x - 1)^4))*cos(2*log(exp(H*(x - 1)^4) - 1) + exp(H*(x - 1)^4))*sqrt(exp(H*(x - 1)^4) - 1)*(x - 1)^3)/(exp(H*(x - 1)^4) + 1)^1.5 - (2*H*exp((-(0.5*H*(x - 1)^4)))*cos(2*log(exp(H*(x - 1)^4) - 1) + exp(H*(x - 1)^4))*(x - 1)^3)/((exp(H*(x - 1)^4) - 1)^0.5*(exp(H*(x - 1)^4) + 1)^0.5)) + 4*H*(x - 1)^3*((16*H^2*exp(1.5*H*(x - 1)^4)*cos(2*log(exp(H*(x - 1)^4) - 1) + exp(H*(x - 1)^4))*(exp(H*(x - 1)^4) + 1)^1.5*(x - 1)^6)/(exp(H*(x - 1)^4) - 1)^1.5 + (12*H*exp((0.5*H*(x - 1)^4))*sin(2*log(exp(H*(x - 1)^4) - 1) + exp(H*(x - 1)^4))*sqrt(exp(H*(x - 1)^4) + 1)*(x - 1)^2)/(exp(H*(x - 1)^4) - 1)^0.5 + (6*H*exp((0.5*H*(x - 1)^4))*cos(2*log(exp(H*(x - 1)^4) - 1) + exp(H*(x - 1)^4))*sqrt(exp(H*(x - 1)^4) - 1)*(x - 1)^2)/(exp(H*(x - 1)^4) + 1)^1.5 - (6*H*exp((-(0.5*H*(x - 1)^4)))*cos(2*log(exp(H*(x - 1)^4) - 1) + exp(H*(x - 1)^4))*(x - 1)^2)/((exp(H*(x - 1)^4) - 1)^0.5*(exp(H*(x - 1)^4) + 1)^0.5) + (8*H^2*exp((0.5*H*(x - 1)^4))*sin(2*log(exp(H*(x - 1)^4) - 1) + exp(H*(x - 1)^4))*sqrt(exp(H*(x - 1)^4) + 1)*(x - 1)^6)/(exp(H*(x - 1)^4) - 1)^0.5 + (4*H^2*exp((0.5*H*(x - 1)^4))*cos(2*log(exp(H*(x - 1)^4) - 1) + exp(H*(x - 1)^4))*sqrt(exp(H*(x - 1)^4) - 1)*(x - 1)^6)/(exp(H*(x - 1)^4) + 1)^1.5 + (8*H^2*exp((0.5*H*(x - 1)^4))*sin(2*log(exp(H*(x - 1)^4) - 1) + exp(H*(x - 1)^4))*sqrt(exp(H*(x - 1)^4) + 1)*(x - 1)^6)/(exp(H*(x - 1)^4) - 1)^1.5 - (12*H^2*exp(1.5*H*(x - 1)^4)*cos(2*log(exp(H*(x - 1)^4) - 1) + exp(H*(x - 1)^4))*sqrt(exp(H*(x - 1)^4) - 1)*(x - 1)^6)/(exp(H*(x - 1)^4) + 1)^2.5 + (4*H^2*exp((-(0.5*H*(x - 1)^4)))*cos(2*log(exp(H*(x - 1)^4) - 1) + exp(H*(x - 1)^4))*(x - 1)^6)/((exp(H*(x - 1)^4) - 1)^0.5*(exp(H*(x - 1)^4) + 1)^0.5) - (8*H^2*exp(1.5*H*(x - 1)^4)*sin(2*log(exp(H*(x - 1)^4) - 1) + exp(H*(x - 1)^4))*sqrt(exp(H*(x - 1)^4) + 1)*(x - 1)^6)/(exp(H*(x - 1)^4) - 1)^1.5 + (4*H^2*exp((0.5*H*(x - 1)^4))*cos(2*log(exp(H*(x - 1)^4) - 1) + exp(H*(x - 1)^4))*(x - 1)^6)/((exp(H*(x - 1)^4) - 1)^0.5*(exp(H*(x - 1)^4) + 1)^1.5) + (4*H^2*exp((0.5*H*(x - 1)^4))*cos(2*log(exp(H*(x - 1)^4) - 1) + exp(H*(x - 1)^4))*(x - 1)^6)/((exp(H*(x - 1)^4) - 1)^1.5*(exp(H*(x - 1)^4) + 1)^0.5) + (4*H^2*exp(1.5*H*(x - 1)^4)*cos(2*log(exp(H*(x - 1)^4) - 1) + exp(H*(x - 1)^4))*(x - 1)^6)/((exp(H*(x - 1)^4) - 1)^0.5*(exp(H*(x - 1)^4) + 1)^1.5)) + (12*H*exp((-(0.5*H*(x - 1)^4)))*cos(2*log(exp(H*(x - 1)^4) - 1) + exp(H*(x - 1)^4))*(2*x - 2)*sqrt(exp(H*(x - 1)^4) - 1))/(exp(H*(x - 1)^4) + 1)^0.5 - (24*H^2*exp((-(0.5*H*(x - 1)^4)))*cos(2*log(exp(H*(x - 1)^4) - 1) + exp(H*(x - 1)^4))*sqrt(exp(H*(x - 1)^4) - 1)*(x - 1)^5)/(exp(H*(x - 1)^4) + 1)^0.5 - (48*H^2*exp((0.5*H*(x - 1)^4))*sin(2*log(exp(H*(x - 1)^4) - 1) + exp(H*(x - 1)^4))*sqrt(exp(H*(x - 1)^4) + 1)*(x - 1)^5)/(exp(H*(x - 1)^4) - 1)^0.5 - (24*H^2*exp((0.5*H*(x - 1)^4))*cos(2*log(exp(H*(x - 1)^4) - 1) + exp(H*(x - 1)^4))*sqrt(exp(H*(x - 1)^4) - 1)*(x - 1)^5)/(exp(H*(x - 1)^4) + 1)^1.5 + (24*H^2*exp((0.5*H*(x - 1)^4))*cos(2*log(exp(H*(x - 1)^4) - 1) + exp(H*(x - 1)^4))*(x - 1)^5)/((exp(H*(x - 1)^4) - 1)^0.5*(exp(H*(x - 1)^4) + 1)^0.5)))/(exp(H*(x - 1)^4) + 1)^0.5;
f2f1f1f1 = - f1f2f1f1;
f2f1f1f2 = - f1f2f1f2;
f2f1f2f1 = - f1f2f2f1;
f2f1f2f2 = - f1f2f2f2;

f1f2f1f3 = f1f2f1f1;
f1f2f1f4 = f1f2f1f2;

% f1f2f2f1 = 
% f1f2f2f2 = 
f1f2f2f3 = f1f2f2f1;
f1f2f2f4 = f1f2f2f2;

f1f2f3f1 = f1f2f1f1;
f1f2f3f2 = f1f2f1f2;
f1f2f3f3 = f1f2f1f2;
f1f2f3f4 = f1f2f1f2;

f1f2f4f1 = f1f2f2f1;
f1f2f4f2 = f1f2f2f2;
f1f2f4f3 = f1f2f2f1;
f1f2f4f4 = f1f2f2f2;


% 23
f2f3f1f1 = f2f1f1f1;
f2f3f1f2 = f2f1f1f2;
f2f3f1f3 = f2f1f1f1;
f2f3f1f4 = f2f1f1f2;

f2f3f2f1 = f2f1f2f1;
f2f3f2f2 = f2f1f2f2;
f2f3f2f3 = f2f1f2f1;
f2f3f2f4 = f2f1f2f2;

f2f3f3f1 = f2f1f1f1;
f2f3f3f2 = f2f1f1f2;
f2f3f3f3 = f2f1f1f1;
f2f3f3f4 = f2f1f1f2;

f2f3f4f1 = f2f1f2f1;
f2f3f4f2 = f2f1f2f2;
f2f3f4f3 = f2f1f2f1;
f2f3f4f4 = f2f1f2f2;

%14
f1f4f1f1 = f1f2f1f1;
f1f4f1f2 = f1f2f1f2;
f1f4f1f3 = f1f2f1f1;
f1f4f1f4 = f1f2f1f2;

f1f4f2f1 = f1f2f2f1;
f1f4f2f2 = f1f2f2f2;
f1f4f2f3 = f1f2f2f1;
f1f4f2f4 = f1f2f2f2;

f1f4f3f1 = f1f2f1f1;
f1f4f3f2 = f1f2f1f2;
f1f4f3f3 = f1f2f1f1;
f1f4f3f4 = f1f2f1f2;

f1f4f4f1 = f1f2f2f1;
f1f4f4f2 = f1f2f2f2;
f1f4f4f3 = f1f2f2f1;
f1f4f4f4 = f1f2f2f2;

%34
f3f4f1f1 = f1f2f1f1;
f3f4f1f2 = f1f2f1f2;
f3f4f1f3 = f1f2f1f1;
f3f4f1f4 = f1f2f1f2;

f3f4f2f1 = f1f2f2f1;
f3f4f2f2 = f1f2f2f2;
f3f4f2f3 = f1f2f2f1;
f3f4f2f4 = f1f2f2f2;

f3f4f3f1 = f1f2f1f1;
f3f4f3f2 = f1f2f1f2;
f3f4f3f3 = f1f2f1f1;
f3f4f3f4 = f1f2f1f2;

f3f4f4f1 = f1f2f2f1;
f3f4f4f2 = f1f2f2f2;
f3f4f4f3 = f1f2f2f1;
f3f4f4f4 = f1f2f2f2;


all_vec_field = [f1f2f1f1 ;f1f2f1f2 ;f1f2f1f3 ;f1f2f1f4 ;
f1f2f2f1 ;f1f2f2f2;f1f2f2f3;f1f2f2f4 ;
f1f2f3f1 ;f1f2f3f2 ;f1f2f3f3 ;f1f2f3f4 ;
f1f2f4f1 ;f1f2f4f2 ;f1f2f4f3 ;f1f2f4f4 ;
%
f1f4f1f1 ;f1f4f1f2 ;f1f4f1f3 ;f1f4f1f4 ;
f1f4f2f1 ;f1f4f2f2 ;f1f4f2f3 ;f1f4f2f4 ;
f1f4f3f1 ;f1f4f3f2 ;f1f4f3f3 ;f1f4f3f4 ;
f1f4f4f1 ;f1f4f4f2 ;f1f4f4f3 ;f1f4f4f4 ;
%
f2f3f1f1 ;f2f3f1f2 ;f2f3f1f3 ;f2f3f1f4 ;
f2f3f2f1 ;f2f3f2f2 ;f2f3f2f3 ;f2f3f2f4 ;
f2f3f3f1 ;f2f3f3f2 ;f2f3f3f3 ;f2f3f3f4 ;
f2f3f4f1 ;f2f3f4f2 ;f2f3f4f3 ;f2f3f4f4 ;
%
f3f4f1f1 ;f3f4f1f2 ;f3f4f1f3 ;f3f4f1f4 ;
f3f4f2f1 ;f3f4f2f2 ;f3f4f2f3 ;f3f4f2f4 ;
f3f4f3f1 ;f3f4f3f2 ;f3f4f3f3 ;f3f4f3f4 ;
f3f4f4f1 ;f3f4f4f2 ;f3f4f4f3 ;f3f4f4f4 ];

xdot1 = dynamicsEg4_LBS(x,P);
% xdot1and2 = dynamicsEg4_LBS2ndOrder(x,P);
xdot3 =  sum(all_beta.*all_vec_field);


x_dot = xdot1+xdot3;

out = x_dot;
end
