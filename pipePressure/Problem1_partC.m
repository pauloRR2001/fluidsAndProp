clc
clear 
close all

%% initialization

R1 = 1.5; % m
R2 = 2; % m

w1 = 2000; % rev/min
w2 = 500; % rev/min

P1 = 0; % kPa
rho = 1000; % kg/m^3

r = linspace(R1,R2,1000); % m
theta = linspace(0,2*pi,1000); % rad

%% calculations

% convert to rad/s
w1 = w1*(2*pi)/60; % rad/s
w2 = w2*(2*pi)/60; % rad/s

% polynomial coefficients
g1 =  (w2*R2^2-w1*R1^2)/(R2^2-R1^2);
g2 = ((R1*R2)^2)*(w1-w2)/(R2^2-R1^2);

% theta velocity
v_theta = g1*r+g2./r; % m/s

% P(r)
P = rho*(g1^2)*(r.^2-R1^2)/2 + 2*rho*g1*g2*log(r/R1) - rho*(g2^2)*((1./r.^2)-1/R1^2) + P1; % N/m^2

% tangential velocity of each cylinder
v_pipe_in = R1*w1;
v_pipe_out = R2*w2;

%% plots

figure(1)
plot(r,v_theta)
hold on
plot(R1,v_pipe_in,'r*')
plot(R2,v_pipe_out,'b*')
grid on
title('Theta Velocity vs Radius - Part C')
xlabel('radius (m)')
ylabel('v {theta} (m/s)')
legend('v {theta}','v inner pipe','v outer pipe')

figure(2)
plot(r,P/1000)
grid on
title('Pressure vs Radius - Part C')
xlabel('radius (m)')
ylabel('P(r) (kPa)')

%% extra velocity (mesh format)

r = linspace(R1,R2,10);
theta = linspace(0,2*pi,30);

[r, theta] = meshgrid(r,theta);
x = r.*cos(theta);
y = r.*sin(theta);

v_theta = (w2*R2^2-w1*R1^2)*r/(R2^2-R1^2)+((R1*R2)^2)*(w1-w2)*(1./r)/(R2^2-R1^2); % m/s

u = -sin(theta).*v_theta; % x component of the vector field
v = cos(theta).*v_theta; % y component of the vector field

figure(3)
quiver(x, y, u, v)
title('Velocity Field in Pipe Cross-section (m/s) - Part C')
hold on
rectangle('Position', [-R1, -R1, 2*R1, 2*R1], 'Curvature', [1, 1], 'EdgeColor', 'r');
rectangle('Position', [-R2, -R2, 2*R2, 2*R2], 'Curvature', [1, 1], 'EdgeColor', 'r');
axis equal
grid on
xlabel('x (m)')
ylabel('y (m)')

%% extra pressure (mesh format)

P = rho*(g1^2)*(r.^2-R1^2)/2 + 2*rho*g1*g2*log(r/R1) - rho*(g2^2)*((1./r.^2)-1/R1^2) + P1;

u = -sin(theta).*v_theta; % x component of the vector field
v = cos(theta).*v_theta; % y component of the vector field

figure(4);
pcolor(x, y, P/1000);
title('Pressure Gradient in the Pipe (kPa) - Part C')
shading interp;  % Interpolate colors for smooth representation
colorbar;       % Display color bar
grid on
axis equal
xlabel('x (m)')
ylabel('y (m)')


