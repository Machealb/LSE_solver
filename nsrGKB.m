function [U, Z, B, bbeta] = nsrGKB(A, C, b, k, tol, reorth)
    % Null Space Restricted Golub-Kahan bidiagonalizition,
    % where the right Lanczos vectors are in N(C).
    % 
    % Inputs:
    %   A, C: either (a) a full or sparse mxn matrix;
    %             (b) a matrix object that performs the matrix*vector operation
    %   b: right-hand side vector
    %   k: the maximum number of iterations 
    %   tol: stopping tolerance of pcg.m for solving C*sb = s or min||Csb-s||_2
    %       if tol=0, then solve it directly 
    %   reorth: 
    %       0: no reorthogonalization
    %       1: full reorthogonaliation, MGS
    %       2: double reorthogonaliation, MGS
    %
    % Outputs:
    %   U: mx(k+1) column 2-orthornormal matrix
    %   Z: nx(k+1) matrix, column 2-orthonormal 
    %   bbeta: 2-norm of b
    %   B: (k+1)x(k+1) lower bidiagonal matrix
    %
    % Reference: [1]. Haibo Li, Krylov iterative methods for linear least squares
    % problems with linear equality constraints
    %
    % Haibo Li, School of Mathematics and Statistics, The University of Melbourne
    % 31, Dec, 2024.
    
    % Check for acceptable number of input arguments
    if nargin < 6
        error('Not Enough Inputs')
    end

    if isa(C, 'function_handle') && tol == 0
        error('Tol must not be 0 for a funtional handel C')
    end

    [m, n] = size(A); 
    if size(b,1) ~= m || size(C,2) ~= n
        error('The dimensions are not consistent')
    end

    if tol==0
        P = eye(n) - pinv(full(C))*C;
    end

    % declares the matrix size
    fprintf('[Start nsr-GKB...], max_Iter=%d, reorth=%d\n', [k,reorth]);
    B  = zeros(k+1, k+1);
    U  = zeros(m, k+1);
    Z  = zeros(n, k+1); 

    % initial step of nsr-GKB
    bbeta = norm(b);
    u = b / bbeta;  
    U(:,1) = u;
    
    rb = A' * u; 
    if tol == 0
        r = P * rb;
    else
        t = lsqr(C, C*rb,tol, 2*n);
        r = rb - t;
    end

    alpha = sqrt(r'*r);
    z  = r / alpha;     Z(:,1)  = z;
    B(1,1)  = alpha;

    % k step iteration of nsr-GKB
    for j = 1:k
        fprintf('[nsr-GKB iterating...], step=%d--------\n', j);
        % compute u in 2-inner product
        s = A * z - alpha * u;
        if reorth == 1
            for i = 1:j
                s = s - U(:,i)*(U(:,i)'*s);  % MGS
            end
        elseif reorth == 2
            for i = 1:j
                s = s - U(:,i)*(U(:,i)'*s);  % MGS
            end
            for i = 1:j
                s = s - U(:,i)*(U(:,i)'*s);  % MGS
            end
        end

        beta = norm(s);
        if beta < 1e-14
            fprintf('[Breakdown...], beta=%f, nsr-GKB breakdown at %d--\n', [beta,j]);
            U  = U(:,1:j);
            Z  = Z(:,1:j);
            B  = B(1:j,1:j);
            break;
        end

        u = s / beta;
        U(:,j+1) = u;
        B(j+1,j) = beta;

        % compute z in 2-inner product
        rb = A' * u;
        if tol == 0
            r = P * rb;
        else
            t = lsqr(C, C*rb,tol, 2*n);
            r = rb - t;
        end

        r = r - beta * z;

        if reorth == 1
            for i = 1:j
                r = r - Z(:,i)*(Z(:,i)'*r);
            end
        elseif reorth == 2
            for i = 1:j
                r = r - Z(:,i)*(Z(:,i)'*r);
            end
            for i = 1:j
                r = r - Z(:,i)*(Z(:,i)'*r);
            end
        end
        
        alpha = sqrt(r'*r);
        if alpha < 1e-14
            fprintf('[Breakdown...], alpha=%f, nsr-GKB breakdown at %d--\n', [alpha,j]);
            U  = U(:,1:j+1);
            Z  = Z(:,1:j);
            B  = B(1:j+1,1:j);
            break;
        end

        z = r / alpha;
        Z(:,j+1) = z;
        B(j+1,j+1) = alpha;
    end

end