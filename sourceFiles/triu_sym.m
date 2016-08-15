function [W] = triu_sym(a)
%TRIUL - Triangular matrix Upper Left (Symbolic).
%
% a         : eigenvector of components
% W         : triangular matrix 
% Author    : Thomas Beauduin, University of Tokyo, 2016
%%%%%
N = length(a);

for r=1:N
    for c=1:N
        if     r+c <= N,   W(r,c) = a(r+c); 
        elseif r+c == N+1, W(r,c) = 1.0;
        elseif r+c >  N+1, W(r,c) = 0.0; 
        end
    end
end

end

