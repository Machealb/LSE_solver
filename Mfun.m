function y = Mfun(z, A, L)
% Let G = A'A+L'*L, 
% This function computes G*z for a vector z.
%
% Haibo Li, Institute of Computing Technology, Chinese Academy of Sciences
% 04, July, 2023.

if isa(A, 'function_handle')
    az = A(z, 'notransp');
    v1 = A(az, 'transp');
else
    v1 = A' * (A * z);
end

if isa(L, 'function_handle')
    az = L(z, 'notransp');
    v2 = L(az, 'transp');
else
    v2 = L' * (L * z);
end

y = v1 + v2;

end
