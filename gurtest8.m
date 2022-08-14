function [x y]=gurtest8(in)
%If system of equations has no solution, we can use this function
%to find the new modification factors F' for all f. 

%AMIN: 
%The objective:
%min(x1+x2+...+y1+y2+...+(newf1)^2-(oldf1)*(newf1)+(newf2)^2-(oldf2)*(newf2)+...)

[num_eq,~]=size(in);
x_var=unique(in(:,[1 4]),'rows');
y_var=unique(in(:,[2 5]),'rows');
%num_var_f=num_eq;

%all_var=[x_var;y_var];
[num_var_x,~]=size(x_var);
[num_var_y,~]=size(y_var);
m.A = sparse(0,num_var_x+num_var_y);
% m.A(1)=1;
% m.sense = '<';
% m.rhs = 100;%<-------------------------------------------------
m.lb=0.1*ones(1,num_var_x+num_var_y);
% % Objective function 
m.obj = ones(num_var_x+num_var_y,1);
m.obj(1:num_var_x)=-2*in(1:num_var_x,3);
m.obj(num_var_x+1:end)=-2; 
m.Q=sparse(diag(ones(1,num_var_x+num_var_y)));
params.MIPGap = 1e-5;

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

%j=1;
for i=1:num_eq
    var1=in(i,[1 4]);
    var2=in(i,[2 5]);
    m.quadcon(i).Qrow=full(mat1(var1(1),var1(2)));
    m.quadcon(i).Qcol=full(mat2(var2(1),var2(2)));
    m.quadcon(i).Qval = 1.0;
    % m.quadcon(i).q = sparse(num_var_x+num_var_y+num_var_f,1);
    m.quadcon(i).sense = '=';
 %   if in(i,7)==1
%         lcol=i;
%         lval=-in(i,3);
       m.quadcon(i).q=sparse(num_var_x+num_var_y, 1);
        %j=j+1;
        m.quadcon(i).rhs = in(i,7); 
%     else
%         lcol=1;
%         lval=0;
%         m.quadcon(i).q=sparse(lcol, 1, lval, num_var_x+num_var_y+num_var_f, 1);
%         m.quadcon(i).rhs = in(i,3)*in(i,7);
   % end 
end
% Solve bilinear model, display solution.  The problem is non-convex,
% we need to set the parameter 'NonConvex' in order to solve it.
params.NonConvex = 2;
result = gurobi(m, params);
%disp(result.x);
%y=result;
%toc
x=result.x(1:num_var_x);
y=result.x(num_var_x+1:end);
%f=result.x(num_var_x+num_var_y+1:end);
% Constrain 'x' to be integral and solve again
% m.vtype = 'ICC';
% result = gurobi(m, params);
% disp(result.x);

