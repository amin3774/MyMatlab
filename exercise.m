syms lamda t
Q=lamda.*[-3 3 0 0;0 -2 2 0;0 0 -1 1;0 0 0 0];
pi0=[1 0 0 0];
pi=pi0*expm(Q*t);
transpose(pi)