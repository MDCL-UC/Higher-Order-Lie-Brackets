%% For States
yyaxis left
fontsize(10,"points")


% For 4 input system
sum_val_4i = cumtrapz(t_x_4i,abs(X_4i(:,1).^2));
plot(t_x_4i,sum_val_4i,"color",[0.3010 0.7450 0.9330],"LineWidth",2)
hold on
sum_state_4i = sum_val_4i(end)
ylabel("State effort")

% For Newton based ESC
load NewtonBased2019_4thordercost
sum_val = cumtrapz(t_Newton_based_ESC,abs(x1_Newton_based_ESC).^2);
plot(t_Newton_based_ESC,sum_val,"color",[0.3010 0.7450 0.9330],"LineStyle","--","LineWidth",2)
sum_state_Newton = sum_val(end)
hold on

%% Control inputs
% For 4 control input system
yyaxis right
t = 0:0.1:P.simruntime;
u1 = P.w^P.p1*sin(P.k1*P.w*t);
u2 = P.w^P.p2*cos(P.k2*P.w*t);
u3 = P.w^P.p3*sin(P.k3*P.w*t);
u4 = P.w^P.p4*cos(P.k4*P.w*t);

sum_val1 = cumtrapz(t,abs(u1).^2);
sum_val2 = cumtrapz(t,abs(u2).^2);
sum_val3 = cumtrapz(t,abs(u3).^2);
sum_val4 = cumtrapz(t,abs(u4).^2);

sum_val = sum_val1 + sum_val2+ sum_val3+ sum_val4;
plot(t,sum_val,'r-',"LineWidth",2)
hold on
sum_u_4i = sum_val(end)


% Newton based system
w = 100; % rad/s
k1 = 1;
k2 = 2;
k = k1;
k_z = 2*k1;
a_z = 8*k^2;
p_star = 0.51;

u1 = w^p_star*sin(k*w*t);
u2 = w^(1-p_star)*cos(k*w*t);
u3 = w^(2-2*p_star)* cos(k_z*w*t);

sum_val1 = cumtrapz(t,abs(u1).^2);
sum_val2 = cumtrapz(t,abs(u2).^2);
sum_val3 = cumtrapz(t,abs(u3).^2);

sum_val = sum_val1 + sum_val2+ sum_val3;
plot(t,sum_val,'r--',"LineWidth",2)
sum_u_Newton =  sum_val(end)
ylabel("Control effort")
legend("State effort for proposed ESC","State effort for ESC in [1]","Control effort for proposed ESC","Control effort for ESC in [1]")
grid on
xlabel("Time")
title("Effort vs time")

fontsize(10,"points")
left_color = [0 0 1];
right_color = [1 0 0 ];
set(fig,'defaultAxesColorOrder',[left_color; right_color]);



