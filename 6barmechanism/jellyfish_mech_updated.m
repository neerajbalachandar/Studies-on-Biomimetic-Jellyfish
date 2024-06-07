clc; clear; close all;
% Constants
N = 100;
l1 = 0.5;
l2 = 0.5;
l3 = 0.5;
w3 = 0.5;
w4 = 0.7;
l4 = 1.1;
l5 = 1;
theta1 = 0.785; % 45 degrees in radians
alpha4 = 1.703; % 98 degrees in radians
alpha3 = 1.555; % 89 degrees in radians

% Slider motion ranges
s1 = linspace(0.2, 0.8, N); % First range
s2 = linspace(0.3, 1.1, N); % Second range

theta3_1 = zeros(1, N);
theta2_1 = zeros(1, N);
theta3_2 = zeros(1, N);
theta2_2 = zeros(1, N);
theta2_and_3 = [0.2897, 1.5533]; % Initial guess for theta2 and theta3

s_offset=0.035;

% Solve for theta2 and theta3 for the first range of slider positions
for i = 1:N
    xsol = fsolve(@(x) loop_closure_1(x, l1, l2, l3, l4, l5, s1(i), theta1), theta2_and_3);
    theta2_1(i) = xsol(1);
    theta3_1(i) = xsol(2);
    theta2_and_3 = [theta2_1(i), theta3_1(i)];
end

% Reset initial guess for the second range
theta2_and_3 = [0.2897, 1.5533];

% Solve for theta2 and theta3 for the second range of slider positions
for i = 1:N
    xsol = fsolve(@(x) loop_closure_2(x, l1, l2, l3, l4, l5, s2(i), theta1), theta2_and_3);
    theta2_2(i) = xsol(1);
    theta3_2(i) = xsol(2);
    theta2_and_3 = [theta2_2(i), theta3_2(i)];
end

% Calculate coupler curve coordinates for both ranges
coupler_curve_x_coord_1 = s1 + l1 * cos(theta1) + l2 * cos(theta2_1) + w3 * cos(theta3_1 - alpha3);
coupler_curve_y_coord_1 = l1 * sin(theta1) + l2 * sin(theta2_1) + w3 * sin(theta3_1 - alpha3);
coupler_curve_x_coord_2 = s2 + l1 * cos(theta1) + l2 * cos(theta2_2) + w3 * cos(theta3_2 - alpha3);
coupler_curve_y_coord_2 = l1 * sin(theta1) + l2 * sin(theta2_2) + w3 * sin(theta3_2 - alpha3);

% Plot the mechanisms for each slider position
for i = 1:N
    % First mechanism coordinates
    x_coord_1 = [s1(i), s1(i) + l1 * cos(theta1), s1(i) + l1 * cos(theta1) + l2 * cos(theta2_1(i)), s1(i) + l1 * cos(theta1) + l2 * cos(theta2_1(i)) + l3 * cos(theta3_1(i)), l4, 0, 0];
    y_coord_1 = [0, l1 * sin(theta1), l1 * sin(theta1) + l2 * sin(theta2_1(i)), l1 * sin(theta1) + l2 * sin(theta2_1(i)) + l3 * sin(theta3_1(i)), l5, l5, 0];
    x_coord_link_4_1 = [s1(i) + l1 * cos(theta1) + l2 * cos(theta2_1(i)), s1(i) + l1 * cos(theta1) + l2 * cos(theta2_1(i)) + w3 * cos(theta3_1(i) - alpha3), l4, s1(i) + l1 * cos(theta1) + l2 * cos(theta2_1(i))];
    y_coord_link_4_1 = [l1 * sin(theta1) + l2 * sin(theta2_1(i)), l1 * sin(theta1) + l2 * sin(theta2_1(i)) + w3 * sin(theta3_1(i) - alpha3), l5, l1 * sin(theta1) + l2 * sin(theta2_1(i))];
    slider_x_1 = [0, s1(i)];
    slider_y_1 = [0, 0];
    box_x_1=[s1(i)-s_offset,s1(i)+s_offset,s1(i)+s_offset,s1(i)-s_offset,s1(i)-s_offset];
    box_y_1=[-s_offset,-s_offset,s_offset,s_offset,-s_offset];

    % Second mechanism coordinates
    x_coord_2 = [s2(i), s2(i) + l1 * cos(theta1), s2(i) + l1 * cos(theta1) + l2 * cos(theta2_2(i)), s2(i) + l1 * cos(theta1) + l2 * cos(theta2_2(i)) + l3 * cos(theta3_2(i)), l4-l5*sind(10), -l5*sind(10), 0];
    y_coord_2 = [0, l1 * sin(theta1), l1 * sin(theta1) + l2 * sin(theta2_2(i)), l1 * sin(theta1) + l2 * sin(theta2_2(i)) + l3 * sin(theta3_2(i)), l5*cosd(10), l5*cosd(10), 0];
    x_coord_link_4_2 = [s2(i) + l1 * cos(theta1) + l2 * cos(theta2_2(i)), s2(i) + l1 * cos(theta1) + l2 * cos(theta2_2(i)) + w3 * cos(theta3_2(i) - alpha3), l4-l5*sind(10), s2(i) + l1 * cos(theta1) + l2 * cos(theta2_2(i))];
    y_coord_link_4_2 = [l1 * sin(theta1) + l2 * sin(theta2_2(i)), l1 * sin(theta1) + l2 * sin(theta2_2(i)) + w3 * sin(theta3_2(i) - alpha3), l5*cosd(10), l1 * sin(theta1) + l2 * sin(theta2_2(i))];
    slider_x_2 = [0, s2(i)];
    slider_y_2 = [0, 0];
    box_x_2=[s2(i)-s_offset,s2(i)+s_offset,s2(i)+s_offset,s2(i)-s_offset,s2(i)-s_offset];
    box_y_2=[-s_offset,-s_offset,s_offset,s_offset,-s_offset];


    % Plot first mechanism
    plot(x_coord_1, y_coord_1, 'k-','LineWidth',1.5)
    hold on;
    plot(x_coord_1, y_coord_1, 'ro')
    plot(x_coord_link_4_1, y_coord_link_4_1, 'k-','LineWidth',1.5)
    plot(slider_x_1, slider_y_1, 'ro')
    plot(coupler_curve_x_coord_1(1:i), coupler_curve_y_coord_1(1:i), 'r--',"LineWidth",3)
    plot(box_x_1,box_y_1,'k-')
    grid off

    % Plot second mechanism
    plot(x_coord_2, y_coord_2, 'b-','LineWidth',1.5)
    plot(x_coord_2, y_coord_2, 'go')
    plot(x_coord_link_4_2, y_coord_link_4_2, 'b-','LineWidth',1.5)
    plot(slider_x_2, slider_y_2, 'go')
    plot(coupler_curve_x_coord_2(1:i), coupler_curve_y_coord_2(1:i), 'g--','LineWidth',3)
    plot(box_x_2,box_y_2,'k-')
    hold off;
    grid off;
    axis equal;
    xlim([-2 * (l2 + l3), 3 * (l2 + l3)]);
    ylim([-2 * (l2 + l3), 3 * (l2 + l3)]);
    title('Jellyfish 6bar mechanism', 50, 'FontName', 'Palatino Linotype')
    xlabel('$x$', 'interpreter', 'latex', 'FontSize', 50, 'FontName', 'Palatino Linotype')
    ylabel('$y$', 'interpreter', 'latex', 'FontSize', 50, 'FontName', 'Palatino Linotype')
    pause(0.01);
end
hold on
plot(1.53255,1.33784,'ko','MarkerSize',10,'MarkerFaceColor','g')
plot(1.79137,1.11877,'ko','MarkerSize',10,'MarkerFaceColor','r')

function F = loop_closure_1(x, l1, l2, l3, l4, l5, s, theta1)
    F(1) = s + l1 * cos(theta1) + l2 * cos(x(1)) + l3 * cos(x(2)) - l4;
    F(2) = l1 * sin(theta1) + l2 * sin(x(1)) + l3 * sin(x(2)) - l5;
end

function F = loop_closure_2(x, l1, l2, l3, l4, l5, s, theta1)
    F(1) = s + l1 * cos(theta1) + l2 * cos(x(1)) + l3 * cos(x(2)) - l4 - l5*cosd(100);
    F(2) = l1 * sin(theta1) + l2 * sin(x(1)) + l3 * sin(x(2)) - l5*sind(100);
end