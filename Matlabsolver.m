% clear
% problem.options = optimoptions('fsolve','Display','iter','algorithm','levenberg-marquardt');
% problem.objective = @myfun;
% problem.x0 =1*ones(7,1);
% problem.solver = 'fsolve';
% x = fsolve(problem)

%%%%%
% function Matlabsolver
clc;clear
%tic 0.497688 seconds
tic
x0 = 10*ones(7,1);
A = [];
b = [];
Aeq = [];
beq = [];
lb = 0.1*ones(7,1);
ub = [];
x = fmincon(@minf,x0,A,b,Aeq,beq,lb,ub,@myfun)
toc
%x(8:17)

function [c,ceq] = myfun(x)
%[num_e , ~]=size(in);
% t1=unique(input(:,[1 4]),'rows');
% t2=unique(input(:,[2 5]),'rows');
%display(in)

ceq = [  x(1) * x(4) - 4*2;
         x(1) * x(5) - 4*3;
         x(1) * x(6) - 4*1;
         x(1) * x(7) - 4*1;

         x(2) * x(4) - 5*1;
         x(2) * x(5) - 5*1.5;
         x(2) * x(6) - 5*.5;
         x(2) * x(7) - 5*.5;

         x(3) * x(4) - 7;
         x(3) * x(5) - 7*1.5;
         x(3) * x(6) - 7*.5;
         x(3) * x(7) - 7*.5];
 c = [];
end
% ceq = [  x(1) * x(4) - 4*x(8);
%          x(1) * x(5) - 4*x(9);
%          x(1) * x(6) - 4*x(10);
%          x(1) * x(7) - 4*x(11);
% 
%          x(2) * x(4) - 5*x(12);
%          x(2) * x(5) - 5*x(13);
%          x(2) * x(6) - 5*5;
%          x(2) * x(7) - 5*3;
% 
%          x(3) * x(4) - 7*x(14);
%          x(3) * x(5) - 7*x(15);
%          x(3) * x(6) - 7*x(16);
%          x(3) * x(7) - 7*x(17)];
%  c = [];
function F = minf(x)
%F = sum(x); % minimise x
F=norm(x);
%F=abs(x(1)-4)+abs(x(2)-5)+abs(x(3)-7)+abs(x(4)-1)+abs(x(5)-1)+abs(x(6)-1)+abs(x(7)-1);
%F = sum(abs(x(8:17)-1)); % we want to have x(8:19) close to 1
end
 
