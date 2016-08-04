function out = sym2tf(tf_sym,Ts)
%sym2tf - Transform tf in symbolic to tf
%
% out = sym2tf(tf_sym,Ts)
% tf_sym: Transfer function in symbolic 
% Ts    : sampling period [s]
% out   : Transfer function
% Author: Wataru Ohnishi, 2016

[N_sym,D_sym] = numden(tf_sym);

N = sym2poly(N_sym);
D = sym2poly(D_sym);

switch nargin
    case 1 % continuous
        out = minreal(tf(N,D));
        
    case 2 % discrete
        out = minreal(tf(N,D,Ts));
end

%{
syms s_sym
syms z_sym

Ts = 100e-6;
sym2tf(1/(s_sym+1))

sym2tf(0.1/(z_sym-0.9),Ts)
%}