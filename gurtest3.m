clc;clear
load indice2
[num_eq,~]=size(in);
x_var=unique(in(:,[1 4]),'rows');
y_var=unique(in(:,[2 5]),'rows');
f_var=find(in(:,7)~=1);
num_var_f=length(f_var);

%all_var=[x_var;y_var];
[num_var_x,~]=size(x_var);
[num_var_y,~]=size(y_var);
% Define data
qrow = [1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3];
qcol = [4, 5, 6, 7, 4, 5, 6, 7, 4, 5, 6, 7];
qval = ones(12);
lcol = [8, 9, 10, 11, 12, 13, 1, 1, 14, 15, 16, 17];
lval = [-4, -4, -4, -4, -5, -5, 0, 0, -7, -7, -7, -7];
rhs = [0, 0, 0, 0, 0, 0, 25, 15, 0, 0, 0, 0];
n = 17;

% Add quadratic constraints
for i=1:length(qrow)
    disp(i);
    m.quadcon(i).Qrow = qrow(i);
    m.quadcon(i).Qcol = qcol(i);
    m.quadcon(i).Qval = qval(i);
    m.quadcon(i).q = sparse(lcol(i), 1, lval(i), n, 1);
    m.quadcon(i).rhs = rhs(i);
    m.quadcon(i).sense = '=';
    m.quadcon(i).name = sprintf('qcon%d', i);
end

% Add variable names
vnames = cell(n,1);
for i=1:n
    vnames{i} = sprintf('x%d', i);
end
m.varnames = vnames;

% No linear constraints
m.A = sparse(0,n);
m.lb=0.1*ones(1,17);
% Solve model and display solution
params.NonConvex = 2;
result = gurobi(m, params);
disp(result.x);