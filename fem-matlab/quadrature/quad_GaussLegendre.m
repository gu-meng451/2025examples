function quad_rule = quad_GaussLegendre(num_GaussPoints, num_Dimensions)

%%
nt = num_GaussPoints;
ned = num_Dimensions;

% use the Golub-Welsch Algorithm
beta = 0.5 ./ sqrt(1 - (2 * (1:nt-1)).^(-2));
T = diag(beta,-1) + diag(beta,1);
[V,D] = eig(T);
xi = diag(D);
w = 2*V(1,:).^2;

%% Form output

quad_rule.method = 'GaussLegendre';
quad_rule.order  = 2*nt-1;

switch ned
    case 1
        quad_rule.nt = nt;
        quad_rule.xi = xi;
        quad_rule.w = w;
        quad_rule.domain = 'line [-1,1]';
        
    case 2
        quad_rule.nt = nt^2;
        quad_rule.w = zeros(quad_rule.nt,1);
        quad_rule.xi = zeros(quad_rule.nt,1);
        quad_rule.eta = zeros(quad_rule.nt,1);
        quad_rule.domain = 'square [-1,1]';
        
        q = 0;
        for i = 1:nt
            for j = 1:nt
                q = q+1;
                quad_rule.w(q) = w(i)*w(j);
                quad_rule.xi(q) = xi(i);
                quad_rule.eta(q) = xi(j);
            end
        end
        
    case 3
        error("3D integration is not implemented yet.")
        
end