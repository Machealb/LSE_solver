function [x1, X2, res1, res2] = KIDS2(A, C, b, d, tol1, tol2, k2)
    % Krylov Iterative Decomposition Solver-II (KIDS-II) for the 
    % Linear least squares problem with linear equaily constraints:
    %    min ||Ax-b||_2   s.t.   ||Cx-d||_2=min .
    %
    % Inputs:
    %   A, C: either (a) a full or sparse mxn matrix;
    %             (b) a matrix object that 
    %   tol1, tol2: stopping tolerance for the inner iteraitonsperforms the matrix*vector operation
    %   b, d: right-hand side vector
    %   k2: the maximum number of iterations for NSR-LSQR
    %
    % Outputs:
    %   x1: solution for min||Cx-d||
    %   X2: store the first k iterative solution for nsr-LSQR
    %   res1:  relative residual norm for x1
    %   res2: relative residual norm for A_{N(C)}^{\dag}\tilde{b}
    %
    % Reference: 
    % [1]. Krylov iterative methods for linear least squares problems with 
    %      linear equality constraints
    % [2]. Haibo Li, A new intepretation of the weighted pseudoinverse and its applications
    %
    % Haibo Li, School of Mathematics and Statistics, The University of Melbourne
    % 31, Dec, 2024.
    
    % Check for acceptable number of input arguments
    if nargin < 7
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

    if tol1 == 0
        x1 = pinv(full(C)) * d;
    else
        x1 = lsqr(C, d,tol1, 1000);
    end
    
    % res1 = norm(C'*(C*x1-d)) / (norm(full(C))*norm(d));
    res1 = 0;
    % b1 = b - A*x1;
    b1 = b - mvp(A,x1);
    [X2, res2, ~] = nsrLSQR(A, C, b1, na, k2, tol2, 1);

end
    
    