function [X, res, bnd] = nsrLSQR(A, C, b, na, k, tol, reorth)
    % Null Space Restricted LSQR algorithm for the generalized LS (GLS) problem:
    %   ||Ax-b||_2=min s.t. x\in N(C)
    %
    % Inputs:
    %   A, C: either (a) a full or sparse mxn matrix;
    %             (b) a matrix object that performs the matrix*vector operation
    %   b: right-hand side vector
    %   na: norm of ||\mathcal{A}||, see [1,Sec 4.2] for details
    %   k: the maximum number of iterations 
    %   tol: stopping tolerance of pcg.m for solving M*sb = s or min||Ms-s||_2
    %       if tol=0, then solve it directly 
    %   reorth: 
    %       0: no reorthogonalization
    %       1: full reorthogonaliation, MGS
    %       2: double reorthogonaliation, MGS
    %
    % Outputs:
    %   X: store the first k iterative solution
    %   res: strore relative residual norms
    %   bnd: upper bound on the residual norms
    %
    % Reference: [1]. Krylov iterative methods for linear least squares 
    %   problems with linear equality constraints
    % [2]. Haibo Li, A new intepretation of the weighted pseudoinverse and its applications
    %
    % Haibo Li, School of Mathematics and Statistics, The University of Melbourne
    % 31, Dec, 2024.
    
    % Check for acceptable number of input arguments
    if nargin < 7
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
    B  = zeros(k+1, k+1);
    U  = zeros(m, k+1);
    Z  = zeros(n, k+1);

    X = zeros(n, k);
    res = zeros(k ,1);  % residual norm
    bnd = zeros(k ,1);  % residual norm

    % initial step of nsr-LSQR
    fprintf('[Start nsr-LSQR...], max_Iter=%d, reorth=%d\n', [k,reorth]);
    bbeta = norm(b);
    Ab = na * bbeta;  % used for computing relative residual norms
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

    % k step iteration of nsr-LSQR
    for j = 1:k 
        fprintf('[nsr-LSQR iterating...], step=%d--------\n', j);
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
        if beta < 1e-16
            fprintf('[Breakdown...], beta=%f, nsr-LSQR breakdown at %d--\n', [beta,j]);
            % U  = U(:,1:j);
            % Z  = Z(:,1:j);
            % B  = B(1:j,1:j);
            X  = X(:,1:j);
            res = res(:,1:j);
            bnd = bnd(:,1:j);
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
        if alpha < 1e-16
            fprintf('[Breakdown...], alpha=%f, nsr-LSQR breakdown at %d--\n', [alpha,j]);
            % U  = U(:,1:j+1);
            % Z  = Z(:,1:j);
            % B  = B(1:j+1,1:j);
            X  = X(:,1:j);
            res = res(:,1:j);
            bnd = bnd(:,1:j);
            break;
        end

        z = r / alpha;
        Z(:,j+1) = z;
        B(j+1,j+1) = alpha;

        % Construct and apply orthogonal transformation
        % rrho = sqrt(rho_bar^2 + beta^2); 
        % c1 = rho_bar / rrho;
        % s1 = beta / rrho; 
        % theta = s1 * alpha; 
        % rho_bar = -c1 * alpha;
        % phi = c1 * phi_bar;
        % phi_bar = s1 * phi_bar;  

        % Update the solution and w_i
        % x = x + (phi/rrho) * w;  
        % w = z - (theta/rrho) * w;

        yj = B(1:j+1,1:j) \ (bbeta*[1;zeros(j,1)]);
        x  = Z(:,1:j) * yj;
        X(:,j) = x;

        % Compute the relative residual norm and its upper bound
        rr = A' * (A*x-b);
        if tol == 0
            ss = P * rr;
        else
            t = lsqr(C, C*rr,tol, 2*n);
            ss = rr - t;
        end

        rn = sqrt(ss'*ss);
        res(j) = rn / Ab; 
        bnd(j) = alpha*beta*abs(yj(j)) / Ab;
        
    end
    
end
    
    