clc;clear all
Q=[-30  10  20  0   0 ; 
    20 -40  5   15  0 ; 
    40  5  -65  10  10;
    25  20  0  -55  10;
    30  0   30  0  -60];

p0=[0.1  0.3  0.25  0.25  0.1];
%p0=[0.2  0.2  0.2   0.2   0.2];
dt=0.001;
%p0=stationaryP';
Mpower=Q*dt+eye(5);
R=[  0  10 20  0  0; 
    20   0  5 15  0; 
    40   5  0 10 10;
    25  20  0  0 10;
    30   0 30  0  0];

D=diag((abs(diag(Q))));

Mjacobi=R/D; %R*D^-1
p_power=p0;
p_jacobi=p0;
p_real=p0;
n=20;
for i=1:n
    p_power(i+1,:)=p_power(i,:)*Mpower;
    p_jacobi(i+1,:)=p_jacobi(i,:)*Mjacobi;
    p_real(i+1,:) = p_real(i,:) * expm(Q*i);
end
disp('Power Method: ')
disp(p_power(n,:))
disp('Jacobi Method:')
disp( p_jacobi(n,:))
disp('Real Value: ')
disp(p_real(n,:))

subplot(1,3,1)
plot([0:n],p_power,'DisplayName','p_power','linewidth',2)
xlabel('iteration')
ylabel('Prob.')
title('Power Method')
legend('pi1','pi2','pi3','pi4','pi5')
ylim([0 0.8])
grid on

subplot(1,3,2)
plot([0:n],p_jacobi,'DisplayName','p_jacobi','linewidth',2)
xlabel('iteration')
ylabel('Prob.')
title('Jacobi Method')
legend('pi1','pi2','pi3','pi4','pi5')
ylim([0 0.8])
grid on

subplot(1,3,3)
plot([0:n],p_real,'DisplayName','p_real','linewidth',2)
xlabel('iteration')
ylabel('Prob.')
title('Real Value')
legend('pi1','pi2','pi3','pi4','pi5')
ylim([0 0.8])
grid on