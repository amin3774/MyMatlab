function y=gurtest4(in)
%If system of equations has no solution, we can use this function
%to find the new modification factors F' for f==1. 

%AMIN: 
%The objective:
%min(x1+x2+...+y1+y2+...+f1+f2+...)




% Linear constraint matrix
% clc;clear
% load in
[num_eq,~]=size(in);
x_var=unique(in(:,[1 4]),'rows');
y_var=unique(in(:,[2 5]),'rows');
 f_var=find(in(:,7)==1);%%%%%%%%%%
 %num_var_f=num_eq;%%%%%%%%%%
num_var_f=length(f_var);
%all_var=[x_var;y_var];
[num_var_x,~]=size(x_var);
[num_var_y,~]=size(y_var);
m.A = sparse(0,num_var_x+num_var_y+num_var_f);% no linear constrains
m.lb=1*ones(1,num_var_x+num_var_y+num_var_f);

% % Objective function minimize(x1+x2+...+y1+y2+...)
m.obj = ones(num_var_x+num_var_y+num_var_f,1);
%m.modelsense = 'min';

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
%%

 j=1;
for i=1:num_eq
    var1=in(i,[1 4]);
    var2=in(i,[2 5]);
    m.quadcon(i).Qrow=full(mat1(var1(1),var1(2)));
    m.quadcon(i).Qcol=full(mat2(var2(1),var2(2)));
    m.quadcon(i).Qval = 1.0;
    % m.quadcon(i).q = sparse(num_var_x+num_var_y+num_var_f,1);
    m.quadcon(i).sense = '=';
     if in(i,7)==1
        lcol=num_var_x+num_var_y+j;
        lval=-in(i,3);
        m.quadcon(i).q=sparse(lcol, 1, lval, num_var_x+num_var_y+num_var_f, 1);
        j=j+1;
        m.quadcon(i).rhs = 0; 
     else
        lcol=1;
        lval=0;
        m.quadcon(i).q=sparse(lcol, 1, lval, num_var_x+num_var_y+num_var_f, 1);
        m.quadcon(i).rhs = in(i,3)*in(i,7);
   end 
end
% Solve bilinear model, display solution.  The problem is non-convex,
% we need to set the parameter 'NonConvex' in order to solve it.
params.NonConvex = 2;
result = gurobi(m, params);
%disp(result.x);
%toc
y=result;
% Constrain 'x' to be integral and solve again
% m.vtype = 'ICC';
% result = gurobi(m, params);
% disp(result.x);
%%
% % Bilinear equality constraint: x(1) * x(4) - 4*2
% m.quadcon(1).Qrow = 1;
% m.quadcon(1).Qcol = 4;
% m.quadcon(1).Qval = 1.0;
% m.quadcon(1).q = sparse(7,1);
% m.quadcon(1).rhs = 8.0;
% m.quadcon(1).sense = '=';
% m.quadcon(1).name = 'bilinear1';
% % Bilinear equality constraint: x(1) * x(5) - 4*3
% m.quadcon(2).Qrow = 1;
% m.quadcon(2).Qcol = 5;
% m.quadcon(2).Qval = 1.0;
% m.quadcon(2).q = sparse(7,1);
% m.quadcon(2).rhs = 4.0*3;
% m.quadcon(2).sense = '=';
% m.quadcon(2).name = 'bilinear2';
% % Bilinear equality constraint: x(1) * x(6) - 4
% m.quadcon(3).Qrow = 1;
% m.quadcon(3).Qcol = 6;
% m.quadcon(3).Qval = 1.0;
% m.quadcon(3).q = sparse(7,1);
% m.quadcon(3).rhs = 4.0;
% m.quadcon(3).sense = '=';
% m.quadcon(3).name = 'bilinear3';
% % Bilinear equality constraint: x(1) * x(7) - 4
% m.quadcon(4).Qrow = 1;
% m.quadcon(4).Qcol = 7;
% m.quadcon(4).Qval = 1.0;
% m.quadcon(4).q = sparse(7,1);
% m.quadcon(4).rhs = 4.0;
% m.quadcon(4).sense = '=';
% m.quadcon(4).name = 'bilinear4';
% % Bilinear equality constraint: x(2) * x(4) - 5
% m.quadcon(5).Qrow = 2;
% m.quadcon(5).Qcol = 4;
% m.quadcon(5).Qval = 1.0;
% m.quadcon(5).q = sparse(7,1);
% m.quadcon(5).rhs = 5.0;
% m.quadcon(5).sense = '=';
% m.quadcon(5).name = 'bilinear5';
% % Bilinear equality constraint: x(2) * x(5) - 5*1.5
% m.quadcon(6).Qrow = 2;
% m.quadcon(6).Qcol = 5;
% m.quadcon(6).Qval = 1.0;
% m.quadcon(6).q = sparse(7,1);
% m.quadcon(6).rhs = 5.0*1.5;
% m.quadcon(6).sense = '=';
% m.quadcon(6).name = 'bilinear6';
% % Bilinear equality constraint: x(2) * x(6) - 5*.5
% m.quadcon(7).Qrow = 2;
% m.quadcon(7).Qcol = 6;
% m.quadcon(7).Qval = 1.0;
% m.quadcon(7).q = sparse(7,1);
% m.quadcon(7).rhs = 5.0*.5;
% m.quadcon(7).sense = '=';
% m.quadcon(7).name = 'bilinear7';
% % Bilinear equality constraint: x(2) * x(7) - 5*.5
% m.quadcon(8).Qrow = 2;
% m.quadcon(8).Qcol = 7;
% m.quadcon(8).Qval = 1.0;
% m.quadcon(8).q = sparse(7,1);
% m.quadcon(8).rhs = 5.0*.5;
% m.quadcon(8).sense = '=';
% m.quadcon(8).name = 'bilinear8';
% % Bilinear equality constraint: x(3) * x(4) - 7
% m.quadcon(9).Qrow = 3;
% m.quadcon(9).Qcol = 4;
% m.quadcon(9).Qval = 1.0;
% m.quadcon(9).q = sparse(7,1);
% m.quadcon(9).rhs = 7.0;
% m.quadcon(9).sense = '=';
% m.quadcon(9).name = 'bilinear9';
% % Bilinear equality constraint: x(3) * x(5) - 7*1.5
% m.quadcon(10).Qrow = 3;
% m.quadcon(10).Qcol = 5;
% m.quadcon(10).Qval = 1.0;
% m.quadcon(10).q = sparse(7,1);
% m.quadcon(10).rhs = 7.0*1.5;
% m.quadcon(10).sense = '=';
% m.quadcon(10).name = 'bilinear10';
% % Bilinear equality constraint: x(3) * x(6) - 7*0.5
% m.quadcon(11).Qrow = 3;
% m.quadcon(11).Qcol = 6;
% m.quadcon(11).Qval = 1.0;
% m.quadcon(11).q = sparse(7,1);
% m.quadcon(11).rhs = 7.0*0.5;
% m.quadcon(11).sense = '=';
% m.quadcon(11).name = 'bilinear11';
% % Bilinear equality constraint: x(3) * x(7) - 7*.5
% m.quadcon(12).Qrow = 3;
% m.quadcon(12).Qcol = 7;
% m.quadcon(12).Qval = 1.0;
% m.quadcon(12).q = sparse(7,1);
% m.quadcon(12).rhs = 7.0*.5;
% m.quadcon(12).sense = '=';
% m.quadcon(12).name = 'bilinear12';
