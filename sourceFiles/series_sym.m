function [As,Bs,Cs,Ds] = series_sym(A1,B1,C1,D1,A2,B2,C2,D2)
%SERIES_SYM - series connection of state-space (symbolic).
As = [A1,0*A1 ; B2*C1,A2];
Bs = [B1 ; B2*D1];
Cs = [D2*C1 , C2];
Ds = D1*D2;
end