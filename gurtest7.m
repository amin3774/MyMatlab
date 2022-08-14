function y=gurtest7(in,R1orig,R2orig)
%Amin: This function is perpared for the example in section
%"What if the system of equations has no solution?".
%We look for a solution where $x$ and $y$ are
%close to their original values and $z$ is close to 1.
%original rates: ox={4,5,7}, oy={1,1,1,1}
%%%%THIS ONLY WORKS FOR THIS EXAMPLE%%%%%

[num_eq,~]=size(in);
x_var=unique(in(:,[1 4]),'rows');
y_var=unique(in(:,[2 5]),'rows');
f_var=in(:,6);
num_var_f=length(f_var);

%all_var=[x_var;y_var];
[num_var_x,~]=size(x_var);
[num_var_y,~]=size(y_var);
m.A = sparse(0,num_var_x+num_var_y+num_var_f);
% m.A(1)=1;
% m.sense = '<';
% m.rhs = 100;
m.lb=0.1*ones(1,num_var_x+num_var_y+num_var_f);
% % Objective function minimize(x1+x2+...+y1+y2+...)
params.TimeLimit = 30;
% m.obj(1:num_var_x+num_var_y)=-2*in(1:end,3);
for i=1:num_var_x
m.obj(i)=-2*R1orig(x_var(i,1),x_var(i,2));
end
for i=1:num_var_y
m.obj(num_var_x+i)=-2*R2orig(y_var(i,1),y_var(i,2));
end
% m.obj(1)=-2*4;m.obj(2)=-2*5;m.obj(3)=-2*7;
% m.obj(4)=-2*1;m.obj(5)=-2*1;m.obj(6)=-2*1;m.obj(7)=-2*1;
for i=1:num_var_f
m.obj(num_var_x+num_var_y+i)=-2*in(i,6);
end
% m.obj(num_var_x+num_var_y+1:num_var_x+num_var_y+num_var_f)=-2;
m.Q=sparse(diag(ones(1,num_var_x+num_var_y+num_var_f)));


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

% j=1;
for i=1:num_eq
    var1=in(i,[1 4]);
    var2=in(i,[2 5]);
    m.quadcon(i).Qrow=full(mat1(var1(1),var1(2)));
    m.quadcon(i).Qcol=full(mat2(var2(1),var2(2)));
    m.quadcon(i).Qval = 1.0;
    % m.quadcon(i).q = sparse(num_var_x+num_var_y+num_var_f,1);
    m.quadcon(i).sense = '=';
%     if in(i,7)==1
        lcol=num_var_x+num_var_y+i;
        lval=-in(i,3);
        m.quadcon(i).q=sparse(lcol, 1, lval, num_var_x+num_var_y+num_var_f, 1);
%         j=j+1;
        m.quadcon(i).rhs = 0; 
%     else
%         lcol=1;
%         lval=0;
%         m.quadcon(i).q=sparse(lcol, 1, lval, num_var_x+num_var_y+num_var_f, 1);
%         m.quadcon(i).rhs = in(i,3)*in(i,7);
%     end 
end
% Solve bilinear model, display solution.  The problem is non-convex,
% we need to set the parameter 'NonConvex' in order to solve it.
params.NonConvex = 2;
result = gurobi(m, params);
disp(result.x);
%toc
y=result;

