function [A_sym,B_sym,C_sym,D_sym] = ctrlPID_sym(mn, bn, kn, Ktn, fp, Ts)
%ctrlPID_sym - Design PID controller with discretization and realization by Sylvester matrix
%              
%   [A_sym,B_sym,C_sym,D_sym] = ctrlPID_sym(mn, bn, kn, Ktn, fp, Ts)
% Nominal plant:	Pn = Ktn/(mn*s^2 + bn*s + kn) 
% fp		: Desired closed loop pole
% Ts_sym    : Sampling time
% [A_sym,B_sym,C_sym,D_sym] : A,B,C,D metrics
% Ref. G. C. Goodwin, S. F. Graebe, and M. E. Salgado, Control System Design. 2000.
% Author    : Wataru Ohnishi, The University of Tokyo, 2016
%%%%%

Spid = [1     0     0     0     0  ;
        bn/mn 1     0     0     0  ;
        kn/mn bn/mn Ktn/mn  0     0  ;
        0     kn/mn 0     Ktn/mn  0  ;
        0     0     0     0     Ktn/mn];
wp = fp*2*pi;
Apid = [1;4*wp;6*wp^2;4*wp^3;wp^4];
para_pid = Spid\Apid;

ac1 = para_pid(2);
bc2 = para_pid(3);
bc1 = para_pid(4);
bc0 = para_pid(5);

num_c_sym = [bc2 bc1 bc0];
den_c_sym = [1 ac1 0];
Ts_sym = Ts;
[A_sym,B_sym,C_sym,D_sym] = c2d_tustin_sym(num_c_sym,den_c_sym,Ts_sym);


%{
mn = 10;
bn = 2;
kn = 1;
Ktn = 1;
fp = 20;
Ts = 100e-6;

[A,B,C,D] = ctrlPID_sym(mn, bn, kn, Ktn, fp, Ts);
G = ss(double(A),double(B),double(C),double(D),Ts);

Gc = ctrlPID(mn, bn, kn, Ktn, fp);
Gd = c2d(Gc,Ts,'tunstin');
figure; bode(G,Gd);
%}
