function [R] = assemble_rhs(ned, nen, nnp, nel, eltype, ...
    x, IEN, LM, quad_rules, fext)
%%
totaldofs = ned*nnp;
R = zeros(totaldofs,1);

%% loop over all elements
for e = 1:nel

    nen_e = nen(eltype(e));
    nee = ned*nen_e;

    % quad rule for each element type
    el_quad = quad_rules{ eltype(e) };

    % setup xe ye
    a = 1:nen_e;
    A = IEN(a,e);
    Xe = x(A);

    if eltype(e) == 1

        % 2-Node line
        [fe] = el1_fe(Xe, el_quad, fext);

    else
        error('Error: unknown element type\n');
        %break
    end

    % add to global stiffness matrix
    for loop1 = 1:nee
        i = LM(loop1,e);
        R(i) = R(i) + fe(loop1);
    end

end


