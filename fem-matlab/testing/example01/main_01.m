%% main_01
% simple example of Bar1D with constant loading


%% setup the workspace
clear all
close all
clc

%%
% add in Matlab routines from FEM project
ROOTDIR = fullfile('../..');
path(pathdef)
addpath(fullfile(ROOTDIR,'preprocmesh'));
addpath(fullfile(ROOTDIR,'assemble'));
addpath(fullfile(ROOTDIR,'postproc'));
addpath(fullfile(ROOTDIR,'quadrature'));
addpath(fullfile(ROOTDIR,'shapefunctions'));

%%
% std header stuff
ned = 1;
nen = std_element_defs();


%% Define parameters

% Material constants
properties.E = 1;
properties.A = 1;
properties.L = 1;
properties.f0 = 1;

fext = @(x) properties.f0;

%% Make a new FEM mesh
nnp = 5;
nel = nnp-1;
nee = 2;
% all elements are type 1, 2-node lines
eltype(1:nel) = 1;


% connectivity
%  IEN(a,e) = A;
%  A: global node number
%  a: local node number
%  e: element number
IEN = zeros(nee, nel);
for e = 1:nel
   IEN(:,e) = [e, e+1]; 
end
x = linspace(0,properties.L,nnp)';
y = zeros(size(x));
z = zeros(size(x));

%% Start plot
[fig,ax] = init_plot_figure();
ax.XLim = [0, properties.L];
ax.YLim = [-0.1, 1.1];

% plot_node(ax, 1:nnp, x,y,z, 'r');
plot_node_labels(ax, 1:nnp, x,y,z);

plot_element(ax, 1:nel, IEN, eltype, x,y,z, 'surf_alpha', 0.2);
plot_element_labels(ax, IEN, eltype, x, y, z);


%% generate or load gaussian quadrature information
quad_rules = set_integration_rules( eltype );

%% Apply BC's
% generate BC tables for all nodes

% this is flagged if a node has an essential BC
fix = zeros(ned,nnp);
fix(1,1) = 1;

% this is list of the values of that BC
g_list = zeros(ned,nnp);
g_list(1,1) = 0;

%% Construct the IM and LM Matricies
[ID, LM, neq, gg, nee, ng, idx_ff, idx_fr] = build_mesh(...
    nnp, IEN, g_list, nen, ned, nel, eltype, fix);

%% build system matrices
ndofs = neq+ng;
qn = zeros(ndofs,1);

% Build K
[K] = assemble_K(ned, nen, nnp, nel, eltype, ...
    x, IEN, LM, quad_rules, properties);

% Build the right-hand-side
R = assemble_rhs(ned, nen, nnp, nel, eltype, ...
    x, IEN, LM, quad_rules, fext);

%% Apply BC's
totaldofs = ned*nnp;
q = zeros(totaldofs, 1);

idx = find(fix==1);
for i = 1:ng % ng = length(idx)
    
    A = idx(i);
    P = ID(1,A);
    
    q(P) = g_list(A);
    
end

%% Solve for the unknowns
K_ff = K(idx_ff,idx_ff);
K_fr = K(1:neq, idx_fr);

% solve:
q(1:neq) = K_ff\( R(idx_ff) - K_fr*q(idx_fr) );

%% make a new plot to show the solution
[fig2,ax2] = init_plot_figure();
ax2.XLim = [-0.1, properties.L];
ax2.YLim = [-0.1, 0.6];
ax2.YLabel.String = "Displacement u(x)";

u_true = @(x,p) p.f0/(2*p.E*p.A)*(2*p.L - x).*x;
fplot(ax2, @(x) u_true(x,properties), [0,properties.L], ...
    'LineWidth',2)

plot_element_solution(ax2, IEN, ID, eltype, x, q, 'color', 'r');
