function [X1, X2, res1, res2] = KIDS1(A, C, b, d, tol1, tol2, k1, k2, type)
    % Krylov Iterative Decomposition Solver-I (KIDS-I) for the 
    % Linear least squares problem with linear equaily constraints:
    %    min ||Ax-b||_2   s.t.   ||Cx-d||_2=min .
    %
    % Inputs:
    %   A, C: either (a) a full or sparse mxn matrix;
    %             (b) a matrix object that performs the matrix*vector operation
    %   b, d: right-hand side vector
    %   tol1, tol2: stopping tolerance for the inner iteraitons
    %   k1: the maximum number of iterations for gLSQR
    %   k2: the maximum number of iterations for NSR-LSQR
    %   type: 
    %       'posi': M is positive definite  (M=A'*A+C'*C)
    %       'semi': M is positive semidefinite
    %
    % Outputs:
    %   X1: store the first k1 iterative solution for gLSQR
    %   X2: store the first k2 iterative solution for nsr-LSQR
    %   res1: relative residual norm for C_{A}^{\dag}d
    %   res2: relative residual norm for A_{N(C)}^{\dag}b
    %
    % Reference: 
    % [1]. Krylov iterative methods for linear least squares problems with 
    %      linear equality constraints
    % [2]. Haibo Li, A new intepretation of the weighted pseudoinverse and its applications
    %
    % Haibo Li, School of Mathematics and Statistics, The University of Melbourne
    % 31, Dec, 2024.
    
    % Check for acceptable number of input arguments
    if nargin < 9
        error('Not Enough Inputs')
    end

    % [m, n] = size(A); 
    % [p, ~] = size(C);
    % if size(b,1) ~= m || size(C,2) ~= n || size(d,1) ~= p
    %     error('The dimensions are not consistent')
    % end

    % kk = 50;
    % [~, ~, B, ~] = nsrGKB(A, C, ones(m,1), kk, 0, 1);
    % kk1 = size(B,2);
    % na = svd(B(kk1,kk1-1));
    na = 0.01;
    Ca = 1;

    [X1, res1, ~] = gLSQR(C, A, d, Ca, k1, tol1, 1, type);
    [X2, res2, ~] = nsrLSQR(A, C, b, na, k2, tol2, 1);

end
    
    