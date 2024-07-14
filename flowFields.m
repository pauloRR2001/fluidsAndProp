close all
clc
clear

%% constants

A = 1; % 1/s
B = -2; % 1/s

%% vector field

% Define the grid
x = linspace(0.001, 0.999, 25);  % x-values between 0.1 and 1
y = linspace(0.001, 0.999, 25);  % y-values between 0.1 and 1

[X, Y] = meshgrid(x, y);

% Calculate the velocity components u and v
u = A*(X.^2) ./ Y;
v = B * X .* log(Y);

magnitude = sqrt(u.^2 + v.^2);

% Create a figure
figure;

% Plot the vector field using quiver
customScaleFactor = 1;  
cyan = [0.2 0.8 0.8];
h = quiver(X, Y, 0.1*u./magnitude, 0.1*v./magnitude, customScaleFactor, Color=cyan);
set(h, 'AutoScale', 'off');  % Turn off autoscaling 

xlabel('X-axis');
ylabel('Y-axis');
title('Vector Field: u = x^2/y, v = -2x*ln(y)');

% Adjust axis limits for a better view
axis([0, 1, 0, 1]);
grid on;

%% streamline

tp = [0 1 2.5]; % s

T = [0 5]; % s

x0 = 0.1; % m
y0 = 0.4; % m
L = 0.1; % m

hold on

[x_ll, y_ll] = streamline(x0, y0, A, B, tp, T); % lower left
[x_lr, y_lr] = streamline(x0+L, y0, A, B, tp, T); % lower right
[x_ur, y_ur] = streamline(x0+L, y0+L, A, B, tp, T); % upper right
[x_ul, y_ul] = streamline(x0, y0+L, A, B, tp, T); % upper left

axis equal
grid on

%% plot

axis equal
xlim([-0.1 1.1])
ylim([0 1])

for k=1:length(tp)
    plotFourSidedPolygon([x_ll(k) x_lr(k) x_ur(k) x_ul(k) x_ll(k)], [y_ll(k) y_lr(k) y_ur(k) y_ul(k) y_ll(k)])
    calculatePolygonArea([x_ll(k) x_lr(k) x_ur(k) x_ul(k) x_ll(k)], [y_ll(k) y_lr(k) y_ur(k) y_ul(k) y_ll(k)])

    text(x_ll(k), y_ll(k), sprintf('(%0.2f, %0.2f)', x_ll(k), y_ll(k)), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
    text(x_ll(k), y_ll(k), '*','Color','k', 'FontSize', 16);
end

legend('flow field','streamline','','','','element t = 0','element t = 1','element t = 2.5')

%% functions

function [x, y] = streamline(x0, y0, A, B, tp, T)

    f_bar0 = [x0; y0]; % [m m; m/s m/s]
    opt = odeset('RelTol',1e-5,'AbsTol',1e-30);
    [t, f_bar] = ode45(@(t, f_bar) dfdt_bar(t, f_bar,A,B), T, f_bar0, opt);

    plot(f_bar(:,1),f_bar(:,2),'r--',LineWidth=2)

    x = [0 0 0];
    y = [0 0 0];
    for k=1:length(tp)
        [~, index] = min(abs(t - tp(k)));
        x(k) = f_bar(index,1);
        y(k) = f_bar(index,2);
    end

end

function [dfdt_bar] = dfdt_bar(t, f_bar, A, B)

    dfdt_bar = [A*(f_bar(1)^2)/f_bar(2); B*f_bar(1)*log(f_bar(2))];

end

function plotFourSidedPolygon(x, y)
    
    % Create a figure and plot the polygon
    plot(x, y, LineWidth=3); % 'b' specifies a blue line

end

function area = calculatePolygonArea(x, y)

    % Calculate the area using the shoelace formula
    area = 0.5 * abs(x(1)*y(2) + x(2)*y(3) + x(3)*y(4) + x(4)*y(1) ...
                    - x(2)*y(1) - x(3)*y(2) - x(4)*y(3) - x(1)*y(4));
end

