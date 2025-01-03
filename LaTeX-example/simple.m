%% Spinning Bead
% by T. Fitzgerald
% This integrates the model of bead sliding on a hoop spinning at constant
% angular velocity

clear all; 
close all;
clc

%% Parameters
gR = 20;
Omega = 80*(1/60)*(2*pi/1);
tf = 10;

%% Define the equation of motion
f = @(t,y) [y(2); -sin(y(1))*(gR - Omega^2*cos(y(1))) ];

%%
% Initial conditions
y0 = [pi/2; 0];

%%
% Integrate the solution
sol = ode23t( f, [0,tf], y0 );

%% Make a figure
figure
time = linspace(0,tf,1000);
y = deval(sol, time, 1);
plot(time, y*180/pi, 'DisplayName', 'Nonlinear model', 'LineWidth', 2)
xlabel('Time t [s]')
ylabel('Angular position theta [deg]')