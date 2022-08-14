clc;clear
% syms U_i
U_i=1;
% syms t_on
t_on=1;
% syms U_o
U_o=1;
% syms t_poff
t_poff=1;
 syms I_o
% syms I_Lmax
I_Lmax=1;
% syms Ts
Ts=1;
% syms fs
fs=1;
% syms L  % inductance in Henry
L=1;
 syms I_N
 syms D % duty cycle
 syms M % ratio between output voltage and input voltage

% D = t_on/Ts;
% fs = 1/Ts;
% I_N = U_i/(2*fs*L);
% I_Lmax = U_i * t_on / L;
% I_o = I_Lmax*t_poff / (2*Ts);
% t_poff = U_i * t_on / (U_o -U_i);
% M = U_o / U_i;

Eqns=[D == t_on/Ts, fs ==1/Ts, I_N == U_i/(2*fs*L),I_Lmax==U_i * t_on / L,I_o== I_Lmax*t_poff / (2*Ts),t_poff==U_i * t_on / (U_o -U_i),M ==U_o / U_i];
y=solve(Eqns,D,fs,I_N,I_Lmax,I_o,t_poff,M)
% y=solve(Eqns,D,I_N,I_o)

