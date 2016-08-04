function [A_sym,B_sym,C_sym,D_sym] = ctrlLPF2_sym(freqHz,dp,Ts)
%ctrlLPF2_sym - Design 2nd order high pass filter with discretization and realization
%
%   [A_sym,B_sym,C_sym,D_sym] = ctrlLPF2_sym(freqHz,dp,Ts)
% freqHz    : Cut off frequency [Hz]
% dp        : Damping
% Ts_sym    : Sampling time
% [A_sym,B_sym,C_sym,D_sym] : A,B,C,D metrics
% Author    : Wataru Ohnishi, The University of Tokyo, 2016
%%%%%

num_c_sym = [0 0 (freqHz*2*pi)^2];
den_c_sym = [1 2*dp*freqHz*2*pi (freqHz*2*pi)^2];
Ts_sym = Ts;
[A_sym,B_sym,C_sym,D_sym] = c2d_tustin_sym(num_c_sym,den_c_sym,Ts_sym);

%{
freqHz = 100;
Ts = 100e-6;
dp = 1;

[A,B,C,D] = ctrlLPF2_sym(freqHz,dp,Ts);
G = ss(double(A),double(B),double(C),double(D),Ts);
figure; bode(G);
%}
