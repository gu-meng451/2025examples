function hlist = plot_element_solution(axes_handle, IEN, ID, eltype, x, qn, varargin)


%% Required inputs
p = inputParser;
addRequired(p, 'axes_handle', @simple_axes_check);
addRequired(p, 'IEN', @isnumeric);
addRequired(p, 'ID', @isnumeric);
addRequired(p, 'eltype', @isnumeric);
addRequired(p, 'x', @isnumeric);
% addRequired(p, 'y', @isnumeric);
% addRequired(p, 'z', @isnumeric);
addRequired(p, 'qn', @isnumeric);


%% Optional Inputs
addParameter(p, 'e_in',     1:size(IEN,2), @isnumeric);
addParameter(p, 'color',    [0, 0.1, 0.6] );
addParameter(p, 'linestyle', '--');


%% Parse the inputs
parse(p, axes_handle, IEN, ID, eltype, x, qn, varargin{:});

e_in      = p.Results.e_in;
linecolor = p.Results.color;
linestyle = p.Results.linestyle;

ned = size(ID,1);

%%
%hnum = 0;
hlist = [];

%%
for e = e_in
    
    if eltype(e) == 1
        % 2-node line
        nen_e = 2;
        A = IEN(1:nen_e,e);
        idx = ID(1,A);
        ue  = qn(idx);
        xe = x(A);
        
        h1 = line(axes_handle, ...
            xe, ue, ...
            'linewidth', 2, 'color', linecolor, 'linestyle', linestyle);

        hlist = [hlist, h1];
          
    else
        fprintf(2',['Error: unknown element type:',...
            ' eltype(%d) = %d\n'],e,eltype(e));
    end
    
end
