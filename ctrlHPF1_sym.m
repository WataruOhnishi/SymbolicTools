function [A_sym,B_sym,C_sym,D_sym] = ctrlHPF1_sym(freqHz,Ts)
%ctrlHPF1_sym - Design 1st order high pass filter with discretization and realization
%
%   [A_sym,B_sym,C_sym,D_sym] = ctrlHPF1_sym(freqHz,Ts)
% freqHz    : Cut off frequency [Hz]
% Ts_sym    : Sampling time
% [A_sym,B_sym,C_sym,D_sym] : A,B,C,D metrics
% Author    : Wataru Ohnishi, The University of Tokyo, 2016
%%%%%

num_c_sym = [1 0];
den_c_sym = [1 freqHz*2*pi];
Ts_sym = Ts;
[A_sym,B_sym,C_sym,D_sym] = c2d_tustin_sym(num_c_sym,den_c_sym,Ts_sym);

%{
freqHz = 100;
Ts = 100e-6;

[A,B,C,D] = ctrlHPF1_sym(freqHz,Ts);
G = ss(double(A),double(B),double(C),double(D),Ts);
figure; bode(G);
%}
