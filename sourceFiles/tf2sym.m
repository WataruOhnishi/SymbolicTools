function out = tf2sym(sys,symval)
%tf2sym - Transform tf to symbolic
%
% out = tf2sym(sys,symval)
% sys   : Transfer function
% symval: symbolic variables in tf
% out   : Transfer function in tf
% Author: Wataru Ohnishi, 2017

[N,D] = tfdata(sys,'v');

if nargin < 2
    symval = 's_sym';
end
s_sym = sym(symval);

syms N_sym D_sym
N_sym = 0;
D_sym = 0;
for k = 1:1:length(N)
    N_sym = N_sym + N(length(N)-k+1)*s_sym^(k-1);
end
for k = 1:1:length(D)
    D_sym = D_sym + D(length(D)-k+1)*s_sym^(k-1);
end

out = N_sym/D_sym;
