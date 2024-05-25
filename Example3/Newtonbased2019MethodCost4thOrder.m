function Newtonbased2019MethodCost4thOrder

%%
clc
clear all
close all




%% Parameters

w = 100; % rad/s
k1 = 1;
k2 = 2;
p_star = 0.51;
rho = 0.3;
w_d = 0.5;
w_y = 20;
w_z = 0.5;
C.H = 1/10;

% Parameters ESC
C.number_of_states = 1;


%% Initial condition
x0 = [4];
d0 = [8/rho]; % to match the initial convergence rate
y0 = [0];
z_vech0 = [0];


%% For Newton Method based ESC
simruntime = 40;
[t_Newton_based_ESC,x_Newton_based_ESC] = ode45(@(t,x) dynamics_Newton_ESC(t,x,w,k1,k2,p_star,rho,w_d,w_y,w_z,@cost_function_2018,C), [0,simruntime],[x0;d0;y0;z_vech0]);
x1_Newton_based_ESC = x_Newton_based_ESC(:,1);

% ODE solver for LBS
[t_Newton_based_LBS,x_Newton_based_LBS] = ode45(@(t,x) dynamics_Newton_LBS(t,x,rho,w_d,w_y,w_z,@gradient_function_2018,@Hessian_function_2018,C), [0,simruntime],[x0;d0;y0;z_vech0]);
x1_Newton_based_LBS = x_Newton_based_LBS(:,1);




%% Figure
figure()
% ESC
plot(t_Newton_based_ESC,x1_Newton_based_ESC,'r');
hold on
% LBS
plot(t_Newton_based_LBS,x1_Newton_based_LBS,'b',"LineWidth",2);
legend("x_{esc}","x_{lbs}")
grid on
% yline(1,"LineWidth",2)
xlabel("time")
ylabel("states")
title("Trajectories of ESC and LBS using proposed Newton based algorithm")



save("NewtonBased2019_4thordercost.mat","x1_Newton_based_ESC","x1_Newton_based_LBS","t_Newton_based_ESC","t_Newton_based_LBS")
% save("x_firstorder.mat","x1_Newton_based_ESC","x1_Newton_based_LBS","t_Newton_based_ESC","t_Newton_based_LBS")
end

function out_val = cost_function_2018(x,C)
H = C.H;

% out_val = 2*(x-1)^2; % Second order cost function
out_val = H*(x-1)^4;
end

function out_val = gradient_function_2018(x,C)
x1 = x(1);
H = C.H;
% gradient = 4*x1-4; % Second order gradient
gradient = 4*H*(x-1)^3;

out_val =gradient;
end

function out_val = Hessian_function_2018(x,C)
x1 = x(1);
H = C.H;
% Hessian = 4; % Second order hessian
Hessian = 12*H*(x-1)^2;

out_val =Hessian;
end


function out_val = P_ij(i,j,n)
out_val = (i-1)*n+ 1.5*i-0.5*i^2+(j-i);
end

function out_val = dynamics_Newton_Classic(~,in_val,rho,receive_gradient_function,receive_hessian_function,C)
% Parameters
% % % % w  k1 = 1 k2  p_star rho  w_d  w_y  w_z  w_v  mu 

n = C.number_of_states;

% The states are x, d, y and z
x = in_val(1:n);

% Value of gradient
gradient_val = receive_gradient_function(x,C);
% Value of Hessian
hessian_val = receive_hessian_function(x,C);

% % Dynamics
% x dot
x_dot = -rho\(hessian_val)*gradient_val ;


out_val = x_dot;

end

function out_val = dynamics_Newton_ESC(t,in_val,w,k1,k2,p_star,rho,w_d,w_y,w_z,receive_cost_function,C)
% Parameters
% % % % w  k1 = 1 k2  p_star rho  w_d  w_y  w_z  w_v  mu 

n = C.number_of_states;

% The states are x, d, y and z
x = in_val(1:n);
d = in_val(n+1:2*n);
y = in_val(2*n+1:3*n);
z_vech = in_val(3*n+1:end);

% Parameter k
if length(z_vech)==1 % Size 1x1 matrix
    k=k1;
elseif length(z_vech)==3
    k = [k1;k2];
end

% Get the cost using cost function
cost_val = receive_cost_function(x,C);


% Getting z from z_vech

if length(z_vech)==1 % Size 1x1 matrix
    z=z_vech;

elseif length(z_vech)==3
    z=[z_vech(1) z_vech(2); z_vech(2)  z_vech(3)];

elseif length(z_vech) ==6
    z=[...
        z_vech(1) z_vech(2) z_vech(3);...
        z_vech(2)  z_vech(4) z_vech(5);...
        z_vech(3) z_vech(5) z_vech(6)];
else
    disp("Error")
end
   
% % Dynamics
% x dot
x_dot = rho*d+ w^p_star*sin(k*w*t); % d and k are vectors


% d dot
d_dot = -w_d*(y+z*d); % The multiplication is matrix  multiplication


% y dot
y_dot = -w_y*(y+2*w^(1-p_star)*cost_val*k.*cos(k*w*t));


if length(z_vech)==1 % Size 1x1 matrix
    z_vech = z;
    a_z = 8*k^2;
    % a_z = 0;
    k_z = 2*k;
    z_vech_dot = -w_z*(z_vech - a_z* w^(2-2*p_star)* cost_val* cos(k_z*w*t));
    
elseif length(z_vech)==3 % Size 2x2 matrix

% vech(z) dot
z_vech1 = z_vech(1);
z_vech2 = z_vech(2);
z_vech3 = z_vech(3);

a_z1 = 8*k1^2;  % Corresponding to position (1,1)
a_z2 = 4*k1*k2; % Corresponding to position (1,2)
a_z3 = 8*k2^2;  % Corresponding to position (2,1)

k_z1 = k1+k1; % Corresponding to position (1,1)
k_z2 = k1+k2; % Corresponding to position (1,2)
k_z3 = k2+k2; % Corresponding to position (2,1)

z_vech_dot1 = -w_z*(z_vech1 - a_z1* w^(2-2*p_star)* cost_val* cos(k_z1*w*t));
z_vech_dot2 = -w_z*(z_vech2 - a_z2* w^(2-2*p_star)* cost_val* cos(k_z2*w*t));
z_vech_dot3 = -w_z*(z_vech3 - a_z3* w^(2-2*p_star)* cost_val* cos(k_z3*w*t));
z_vech_dot = [z_vech_dot1;z_vech_dot2;z_vech_dot3];

else 
    disp("error")
end


out_val = [x_dot;d_dot;y_dot;z_vech_dot];

end



function out_val = dynamics_Newton_LBS(t,in_val,rho,w_d,w_y,w_z,receive_gradient_function,receive_hessian_function,C)
% Parameters
% % % % w  k1 = 1 k2  p_star rho  w_d  w_y  w_z  w_v  mu 

% Number of States
n = C.number_of_states;


% The states are x, d, y and z
x = in_val(1:n);
d = in_val(n+1:2*n);
y = in_val(2*n+1:3*n);
z_vech = in_val(3*n+1:end);




% Getting z from z_vech
if length(z_vech)==1 % Size 1x1 matrix
    z=z_vech;

elseif length(z_vech)==3 % 2x2 matrix
    z=[z_vech(1) z_vech(2); z_vech(2)  z_vech(3)];

elseif length(z_vech) ==6 % 3x3 matrix
    z=[...
        z_vech(1) z_vech(2) z_vech(3);...
        z_vech(2)  z_vech(4) z_vech(5);...
        z_vech(3) z_vech(5) z_vech(6)];
else
    disp("Error")
end

% % Dynamics
% x dot
x_dot = rho*d ;
% x_dot = -rho*receive_gradient_function(x);


% d dot
d_dot = -w_d*(y+z*d); % The multiplication is matrix  multiplication
% d_dot = -w_d*(y); % 



% y dot
y_dot = -w_y*(y-receive_gradient_function(x,C));


% vech(z) dot
% z_dot = z-hessian_val; % original
% z_dot = z-12*(x-1)^2; % corresponding to obj = 4(x-1)^3
% z_dot = z; % If no info from hessian
z_dot = (z-receive_hessian_function(x,C));

if length(z_vech)==1 % Size 1x1 matrix
    z_vech_dot = -w_z*z_dot;

elseif length(z_vech)==3 % 2x2 matrix
    z_vech_dot(1,1) = -w_z*z_dot(1,1); 
    z_vech_dot(2,1) = -w_z*z_dot(1,2); 
    z_vech_dot(3,1) = -w_z*z_dot(2,2); 

elseif length(z_vech) ==6 % 3x3 matrix
    z_vech_dot(1,1) = -w_z*z_dot(1,1); 
    z_vech_dot(2,1) = -w_z*z_dot(1,2); 
    z_vech_dot(3,1) = -w_z*z_dot(1,3); 
    z_vech_dot(4,1) = -w_z*z_dot(2,2); 
    z_vech_dot(5,1) = -w_z*z_dot(2,3); 
    z_vech_dot(6,1) = -w_z*z_dot(3,3); 
end

out_val = [x_dot;d_dot;y_dot;z_vech_dot];

end
