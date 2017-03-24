function out = tf_sym(B_sym,A_sym)
%tf_sym - symbolic version of tf
%
% out = tf_sym(B_sym,A_sym)
% B_sym : numerator polynomial
% A_sym : denominator polynomial
% out   : Transfer function in symbolic
% Author: Wataru Ohnishi, 2017

syms s_sym

nb = length(B_sym);
na = length(A_sym);
s_a = sym('a',[1,na]);
for k = 1:1:na
    s_a(k) = s_sym^(k-1);
end
for k = 1:1:nb
    s_b(k) = s_sym^(k-1);
end
s_a = fliplr(s_a);
s_b = fliplr(s_b);


out = B_sym*s_b.' / (A_sym*s_a.');
