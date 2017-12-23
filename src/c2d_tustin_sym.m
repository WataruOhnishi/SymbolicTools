function [A_sym,B_sym,C_sym,D_sym] = c2d_tustin_sym(num_c_sym,den_c_sym,Ts_sym)
%c2d_tustin_sym - Tustin transform in symbolic
%
%   [A_sym,B_sym,C_sym,D_sym] = func_c2d_tustin_sym(num_c_sym,den_c_sym,Ts_sym)
% num_c_sym : numerator in descending order of powers
% den_c_sym : denominator in descending order of powers
% Ts_sym    : Sampling time
% [A_sym,B_sym,C_sym,D_sym] : A,B,C,D metrics
% Author    : Wataru Ohnishi, The University of Tokyo, 2016
%%%%%

syms s_sym z_sym
n = length(den_c_sym) - 1;
S_sym = sym('s',[n+1,1]);
for kk = 1:1:n+1
    S_sym(kk,1) = s_sym^(kk-1);
end
S_sym = flipud(S_sym);

sys_c = num_c_sym*S_sym/(den_c_sym*S_sym);
% Tustin transform
sys_d = subs(sys_c,{s_sym},{2*(z_sym-1)/(Ts_sym*(z_sym+1))});
temp = children(partfrac(sys_d,z_sym)); % strictly proper
% z_symを含むかどうか
T1 = zeros(length(temp),1); % zの多項式
T2 = ones(length(temp),1); % 定数項
for kk = 1:1:length(temp)
    [num,den] = numden(temp(kk));
    if 1 < length(coeffs(num)) || 1 < length(coeffs(den)) % zの多項式の場合
       T1(kk) = 1;
       T2(kk) = 0;
    end
end
sys_d_proper = temp*T1;
sys_d_const = temp*T2;

[num_d_sym,den_d_sym] = numden(simplifyFraction(sys_d_proper,'Expand',true));
% [num_d_sym,den_d_sym] = numden(temp(1));
if length(coeffs(den_d_sym,'z_sym')) == 1
    [num_d_sym,den_d_sym] = numden(simplifyFraction(temp(2),'Expand',true));
    D_sym = simplify(temp(1));
end
num_d_sym =  coeffs(num_d_sym,'z_sym'); % ascending order of powers
den_d_sym =  coeffs(den_d_sym,'z_sym'); % ascending order of powers

% CCF
A_sym = simplify([zeros(n-1,1) eye(n-1,n-1);
    -den_d_sym(1:n)/den_d_sym(n+1);]);
B_sym = simplify([zeros(n-1,1); 1/den_d_sym(n+1)]);
C_sym = simplify(num_d_sym);
D_sym = simplify(sys_d_const);


%{
% check
Ts_sym = 100e-6;
num_c_sym = [0 2];
den_c_sym = [3 4]; 

num_c_sym = [3 1 2];
den_c_sym = [2 3 4]; 
[A_sym,B_sym,C_sym,D_sym] = c2d_tustin_sym(num_c_sym,den_c_sym,100e-6);
minreal(tf(ss(double(A_sym),double(B_sym),double(C_sym),double(D_sym),100e-6)))
minreal(c2d(tf(num_c_sym,den_c_sym),100e-6,'tustin'))
%}
