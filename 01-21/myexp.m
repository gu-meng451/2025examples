clear all
close all
clc

%%

%%
x = -10;

[y1, iter1] = MyExp(x)

abs(exp(x) - y1)

%%
function [my_exp, iter] = MyExp(x)
iter = 0;
maxiter = 1000;
tol = 1e-8;
flag = 0;
my_exp = 1;
newterm = 1;
while flag == 0
    iter = iter + 1;

    % update exp
    newterm = newterm*x/iter;
    my_exp = my_exp + newterm;

    % absolute error:
    % abs_err = abs( exp(1) - my_exp );

    % rel_err
    err = abs( newterm );

    if iter >= maxiter
        flag = -1;
        error("Failed to converge in %d steps", iter);
    elseif err < tol
        flag = 1;
    end
end

end