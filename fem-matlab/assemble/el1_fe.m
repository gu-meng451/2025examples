function fe = el1_fe(xe, el_quad, external_forcing)

ned = 1;
nen = length(xe);
nee = ned * nen;
fe = zeros(nee,1);

% integration loop
for i = 1:el_quad.nt
    
    xi = el_quad.xi(i);
    w = el_quad.w(i);

    % Evaluate the shape function
    [Ne, Nr] = el1_ShapeFunctions(xi);

    % evaluate the external loading at x(ξ)
    x = dot(Ne, xe);
    fext = external_forcing(x);

    % build Jacobian
    detJ = dot(Nr, xe);

    % dV0
    dV0 = detJ * w;

    % integrate fe:
    fe = fe + Ne' * fext * dV0;
end

end