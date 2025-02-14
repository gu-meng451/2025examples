clear all
close all
clc

%%
k1 = 1;
k2 = 2;
k3 = 1;

m1 = 1;
m2 = 1;

K = [k1+k2, -k2; -k2, k2+k3];
M = diag([m1,m2]);

A = [ zeros(2,2), eye(2);
     -M\K, zeros(2,2)];

B = [ zeros(2,1); M\[1;0] ];

u = @(t) 0;

dx = @(t,x) A*x + B*u(t);

x0 = [0.2; -0.2; 0; 0];

sol = ode45(dx, [0,30],x0);

figure
plot(sol.x, sol.y(1:2,:), "LineWidth",3)

%%
% we want the system to have a settle time of < 0.5 s

sx = [-8+2j, -8-2j, -40, -40];
d = poly(sx)

%%
% let's build P
lambda = eig(A)
a = poly(lambda)

%%
% Atilde (controller canonical form)
n = size(A,1);
Atilde = diag(ones(n-1,1),1);
Atilde(n,:) = -fliplr(a(2:end));

Btilde = zeros(n,1);
Btilde(n) = 1;

Cont_X = [B, A*B, A^2*B, A^3*B]

Cont_Z = [Btilde, Atilde*Btilde, Atilde^2*Btilde, Atilde^3*Btilde]

P = Cont_X/Cont_Z

%%
% controller gains in the z-states
Kz = fliplr( d(2:end) - a(2:end) )

Kx = Kz/P

%% check
eig(A - B*Kx)


%%

u = @(t,x) -Kx*x;

dx = @(t,x) A*x + B*u(t,x);

x0 = [0.2; -0.2; 0; 0];

sol = ode45(dx, [0,1],x0);

figure
plot(sol.x, sol.y(1:2,:), "LineWidth",3)


%%
C = [1, 0, 0, 0];
Ox = [ C; 
       C*A;
       C*A^2;
       C*A^3];
rank(Ox)


