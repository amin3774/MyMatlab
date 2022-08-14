function y=gurtest(in)
%tic 0.210670 seconds
%tic
%AMIN: This function solves the system of 12 nonlinear equations, Minimizes
%the sum of {x1,..,x7} and also tries to minimize each of {x1,..,x7}.

% This example formulates and solves the following simple bilinear model:
%  maximize    x
%  subject to  x + y + z <= 10
%              x * y <= 2         (bilinear inequality)
%              x * z + y * z = 1  (bilinear equality)
%              x, y, z non-negative (x integral in second version)

% Copyright 2020, Gurobi Optimization, LLC
%names = {'x1'; 'x2'; 'x3'; 'x4'; 'x5'; 'x6'; 'x7'};
% Linear constraint matrix
% load indice2
[num_eq,~]=size(in);
x_var=unique(in(:,[1 4]),'rows');
y_var=unique(in(:,[2 5]),'rows');
% all_var=[x_var;y_var];
[num_var_x,~]=size(x_var);
[num_var_y,~]=size(y_var);
 m.A = sparse(ones(1,num_var_x+num_var_y));
% m.A(1)=1;
%m.A = sparse(0,num_var_x+num_var_y);
m.sense = '<';
m.rhs = 10000;
m.lb=0.1*ones(1,num_var_x+num_var_y);
% m.OutputFlag =0;

% Variable names
%m.varnames = names;
m1=max(x_var,[],'all');
m2=max(y_var,[],'all');
mat1=sparse(m1,m1);
mat2=sparse(m2,m2);
idx1=sub2ind([m1 m1],x_var(:,1),x_var(:,2));
idx2=sub2ind([m2 m2],y_var(:,1),y_var(:,2));
mat1(idx1)=[1:num_var_x];
mat2(idx2)=[num_var_x+1:num_var_x+num_var_y];
params.MIPGap = 1e-5;
% Objective function minimize(x1+x2+...+y1+y2+...)
%m.obj = ones(num_var_x+num_var_y,1);
%m.modelsense = 'min';
for i=1:num_eq
    var1=in(i,[1 4]);
    var2=in(i,[2 5]);
    m.quadcon(i).Qrow = full(mat1(var1(1),var1(2)));
    m.quadcon(i).Qcol = full(mat2(var2(1),var2(2)));
    m.quadcon(i).Qval = 1.0;
    m.quadcon(i).q = sparse(num_var_x+num_var_y,1);
    m.quadcon(i).rhs = in(i,7);
    m.quadcon(i).sense = '=';
    %m.quadcon(i).name = 'bilinear1';  
end

% Solve bilinear model, display solution.  The problem is non-convex,
% we need to set the parameter 'NonConvex' in order to solve it.
params.NonConvex = 2;
y = gurobi(m, params);
%disp(result.x);
%toc

end