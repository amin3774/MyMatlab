function [d]=Dependence(B)
%if output is 0 matrix B is independent
% if 1 then dependent.
C=rref(B);
m=length(diag(B(:,1)));
n=length(B(1,:));
if n>m
 d=0;
else
s=sum(diag(C));
if n>s
d=0;
else
d=1;
end
end