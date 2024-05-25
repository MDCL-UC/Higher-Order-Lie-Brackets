figure(1)
%% For States
yyaxis left

% For 4 input system
sum_val_4i = cumtrapz(t_esc4,abs(x_esc4).^2);
plot(t_esc4,sum_val_4i,"color",[0.3010 0.7450 0.9330],"LineWidth",2)
hold on
sum_state_4i = sum_val_4i(end)

% For Grushkovskaya
sum_val = cumtrapz(t_esc2grush,abs(x_esc2grush).^2);
plot(t_esc2grush,sum_val,"color",[0.3010 0.7450 0.9330],"LineStyle","--","LineWidth",2)
hold on
ylabel("State effort")
sum_state_grush = sum_val(end)


%% For control inputs
% For 4 control inputs
yyaxis right

t = 0:0.1:P.simruntime;
u1 = P.w^P.p1*cos(P.k1*P.w*t);
u2 = P.w^P.p2*sin(P.k2*P.w*t);
u3 = P.w^P.p3*cos(P.k3*P.w*t);
u4 = P.w^P.p4*sin(P.k4*P.w*t);

sum_val1 = cumtrapz(t,abs(u1).^2);
sum_val2 = cumtrapz(t,abs(u2).^2);
sum_val3 = cumtrapz(t,abs(u3).^2);
sum_val4 = cumtrapz(t,abs(u4).^2);

sum_val = sum_val1 + sum_val2+ sum_val3+ sum_val4;
plot(t,sum_val,'r',"LineWidth",2)
sum_u_4i = sum_val(end)

% For Grushkovskaya
t = 0:0.1:P.simruntime;
u1 = P.w^0.5*cos(P.k1*P.w*t);
u2 = P.w^0.5*sin(P.k2*P.w*t);

sum_val1 = cumtrapz(t,abs(u1).^2);
sum_val2 = cumtrapz(t,abs(u2).^2);

sum_val = sum_val1 + sum_val2;
plot(t,sum_val,'r--',"LineWidth",2)
sum_u_grush = sum_val(end)

ylabel("Control effort")
legend("State effort for proposed ESC","State effort for ESC in [22]",...
    "Control effort for proposed ESC","Control effort for ESC in [22]")
grid on
xlabel("Time")
title("Effort vs time")


