function [Uc] = ctrbs(A,b)
%CTRBS - Compute Controllability Matrix (Symbolic).
%
% A,b       : symbolic A/b matrices
% Uc        : result
% Author    : Thomas Beauduin, University of Tokyo, 2016
%%%%%
N = length(b);

Uc = b;
for n=2:N
    Uc = [Uc, A^(n-1)*b];
end

end

