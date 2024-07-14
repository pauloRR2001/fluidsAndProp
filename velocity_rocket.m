clc
clear

%% values

m_0 = 10; % kg
m_I = 1; % kg

rho_f = 1.2; % kg/m^3
v_j = 340; % m/s

D_e = 7/100; % m
A_e = (pi/4)*D_e^2; % m^2

c_d = 0.2; % 0.2
A_rf = 120/(100^2); % m^2
rho_a = 1.22; % kg/m^3
g = 9.81; % m/s^2

%% integration with drag

T = [0 15]; % s
X_0 = 0; % [v_0 dv/dt_0], [m/s m/s^2]

opt = odeset('AbsTol', 1e-12, 'RelTol', 1e-8);

[t_no_drag, X_no_drag] = ode45(@(t, X_no_drag) dxdt(t,X_no_drag,rho_a,rho_f,v_j,m_0,A_rf,A_e,0,g), T, X_0,opt);

[t, X] = ode45(@(t, X) dxdt(t,X,rho_a,rho_f,v_j,m_0,A_rf,A_e,c_d,g), T, X_0,opt);

%% analytical without drag

t_f = (m_0-m_I)/(rho_f*v_j*A_e) % s
t_analy = linspace(0,t_f,1000); % s
v_analy = v_j*log(m_0./(m_0-rho_f*v_j*A_e*t_analy)) - g*t_analy; % m/s

%% plot

figure(1)
plot(t,X(:,1), 'r', LineWidth=2)
title('Rocket Velocity vs Time - Problem 3')
ylabel('Rocket Velocity (m/s)')
xlabel('Time (s)')
grid on

hold on
plot(t_analy,v_analy, 'g', LineWidth=2)
plot(t_no_drag,X_no_drag, 'b--', LineWidth=2)

plot([t_f, t_f], ylim, 'k--', 'LineWidth', 1)

legend('Numerical with Drag', 'Analytical without Drag', 'Numerical without Drag', 'Model Limit: m_{fuel}=0')
xlim([0 8])
ylim([0 800])

%% functions

function [dxdt] = dxdt(t,X,rho_a,rho_f,v_j,m_0,A_rf,A_e,c_d,g)

dxdt = [rho_f*(v_j^2)*A_e/(m_0-rho_f*v_j*A_e*t) - rho_a*(X.^2)*A_rf*c_d./(2*(m_0-rho_f*v_j*A_e*t)) - g];

end



