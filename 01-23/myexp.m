clear all
close all
clc

%%
x = 10.1;
N = 1000;

tic
for i = 1:N
    [y1, iter1] = MyExp(x);
end
toc
iter1

abs(exp(x) - y1)

tic
for i = 1:N
    [y2, iter2] = MyExp_iter(x);
end
toc
iter2



tic
for i = 1:N
    y3 = exp(x);
end
toc

%%
function [my_exp, iter] = MyExp(X)

x = X;
if X < 0
    x = -X;
end

% x = n + w
n = floor(x);
w = x - n;

% e^x = e^(n+w) = e^n * e^w
e = 2.718281828459045;
en = 1;
for i = 1:n
    en = en*e;
end

[ew, iter] = MyExp_iter(w);

if X < 0
    my_exp = 1/(en*ew);
else
    my_exp = en*ew;
end

end


%%
function [my_exp, iter] = MyExp_iter(x)
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