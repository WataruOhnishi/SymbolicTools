function [Ad_sym, bd_sym, cd_sym] = c2d_zoh_sym(Ac_sym, bc_sym, cc_sym, Ts_sym, n_expm, n_trun)
%c2d_zoh_sym - c2d by zero-order-hold by symbolic math toolbox
%
% [Ad_sym, bd_sym, cd_sym] = c2d_zoh_sym(Ac_sym, bc_sym, cc_sym, Ts_sym, n_order)
% Ac_sym, bc_sym, cc_sym : A, B, C matrices in continuous
% Ts_sym                 : sampling time [s]
% n_expm                 : maximum order for exact expm
% n_trun                 : truncation order for expm
% Author                 : Wataru Ohnishi, University of Tokyo, 2017
%%%%%

n = length(Ac_sym);
if nargin < 4
    n_expm = 3;
end
if nargin < 5
    n_trun = 10;
end

Ec = sym('Ec_', n+1); % def of augmented system
Ec(1:n, 1:n+1) = [Ac_sym bc_sym]*Ts_sym;
Ec(n+1,:) = zeros(1,n+1);

if n > n_expm
    fprintf('Taylor approximation for n = %d\n',n_trun)
    S = zeros(n+1);
    
    for k = 0:n_trun
        S = S +  Ec^k/factorial(k);
    end
    
    Ad_sym = S(1:n, 1:n);
    bd_sym = S(1:n, n+1);
    
else
    [m,n] = size(Ac_sym); %#ok<ASGLU>
    [m,nb] = size(bc_sym); %#ok<ASGLU>
    s = expm([[Ac_sym bc_sym]*Ts_sym; zeros(nb,n+nb)]);
    Ad_sym = s(1:n,1:n);
    bd_sym = s(1:n,n+1:n+nb);
    % see c2d.m
    cd_sym = cc_sym;
end
