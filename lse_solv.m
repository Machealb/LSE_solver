function [x] = lse_solv(A, C, b, d, method)
    % Direct methods for the Linear least squares problem 
    % with linear equaily constraints:
    %    min ||Ax-b||_2   s.t.   ||Cx-d||_2=min .
    %
    % Inputs:
    %   A, C: either (a) a full or sparse mxn matrix;
    %             (b) a matrix object that performs the matrix*vector operation
    %   b, d: right-hand side vector
    %
    %   method: 
    %       '1': null-space method
    %       '2': direct elimination method
    % Outputs:
    %   x: the solution
    %
    % Reference: 
    % [1]. Krylov iterative methods for linear least squares problems with 
    %      linear equality constraints
    %
    % Haibo Li, School of Mathematics and Statistics, The University of Melbourne
    % 31, Dec, 2024.

    % Check for acceptable number of input arguments
    if nargin < 5
        error('Not Enough Inputs')
    end

    [m, n] = size(A); 
    [p, ~] = size(C);
    if size(b,1) ~= m || size(C,2) ~= n || size(d,1) ~= p
        error('The dimensions are not consistent')
    end

    if strcmp(method, '1')
        x0 = pinv(full(C)) * d;
        V  = null(full(C));
        b1 = b - A * x0;
        A1 = A * V;
        y  = pinv(A1) * b1;
        x = x0 + V * y;
    elseif strcmp(method, '2')
        [U, S, V] = svd(full(C));
        r  = rank(full(S));
        C1 = U(:,1:r) * S(1:r,1:r);
        AA = A*V;
        A1 = AA(:,1:r);
        A2 = AA(:,r+1:n);

        y1 = C1 \ d;
        y2 = A2 \ (b-A1*y1);
        x = V * [y1;y2];
    end

end


