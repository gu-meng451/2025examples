clear all
close all
clc

%%
x = 1;

%%

iter = 0;
maxiter = 1000;
tol = 1e-6;
flag = 0;
my_exp = 1;
newterm = 1;
while flag == 0
    iter = iter + 1;

    % update exp
    newterm = newterm*x/iter;
    my_exp = my_exp + newterm;

    % this does not end yet
end

