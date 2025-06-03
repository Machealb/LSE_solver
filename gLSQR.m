function [X, res, bnd] = gLSQR(A, L, b, Ca, k, tol, reorth, type)
    % Generalized LSQR algorithm for the generalized LS (GLS) problem:
    %    min ||Lx||_2   s.t.   ||Ax-b||_2=min
    % where M = A'*A + L'*L may be positive semi-definite.
    %
    % gLSQR compute the unique minimum 2-norm solution of GLS by
    % applying a Krylov subspace method called the generalized GKB.
    %
    % Inputs:
    %   A, L: either (a) a full or sparse mxn matrix;
    %             (b) a matrix object that performs the matrix*vector operation
    %   b: right-hand side vector
    %   Ca: c-component of the largest generalized singular value of{A, L}
    %   k: the maximum number of iterations 
    %   tol: stopping tolerance of pcg.m for solving M*sb = s or min||Ms-s||_2
    %       if tol=0, then solve it directly 
    %   reorth: 
    %       0: no reorthogonalization
    %       1: full reorthogonaliation, MGS
    %       2: double reorthogonaliation, MGS
    %   type:
    %       'posi': M is positive definite
    %       'semi': M is positive semidefinite
    %
    % Outputs:
    %   X: store the first k iterative solution
    %   res: strore relative residual norms
    %   bnd: upper bound on the residual norms
    %
    % Reference: [1]. Haibo Li, Characterizing GSVD by singular value expansion of
    %   linear operators and its computation
    % [2]. Haibo Li, A new intepretation of the weighted pseudoinverse and its applications
    %
    % Haibo Li, School of Mathematics and Statistics, The University of Melbourne
    % 31, Dec, 2024.
    
    % Check for acceptable number of input arguments
    if nargin < 8
        error('Not Enough Inputs')
    end

    if (isa(A, 'function_handle') || isa(L, 'function_handle')) && tol == 0
        error('Tol must not be 0 for funtional handles')
    end

    [m, n] = sizemm(A);   
    [p,~] = sizemm(L);
    % p = sizemm(L,1);
    % if size(b,1) ~= m || size(L,2) ~= n
    %     error('The dimensions are not consistent')
    % end

    if (isa(A, 'function_handle') || isa(L, 'function_handle'))
        M = @(x) Mfun(x,A,L);
    else
        M = A'*A + L'*L;
    end

    % M1 = M + 1e-3*speye(n);
    % R = ichol(M+1e-12*speye(n));  % incomplete Cholesky to construct a preconditioner for Mx=y
    % Mp = pinv(full(M));
    if strcmp(type, 'semi') && tol == 0
        Mp = pinv(full(M));
    end

    % declares the matrix size
    B  = zeros(k+1, k+1);
    U  = zeros(m, k+1);
    Z  = zeros(n, k+1);
    Zb = zeros(n, k+1);

    X = zeros(n, k);
    res = zeros(k ,1);  % residual norm
    bnd = zeros(k ,1);  % residual norm

    % initial step of gGKB
    fprintf('[Start gLSQR...], max_Iter=%d, reorth=%d\n', [k,reorth]);
    bbeta = norm(b);
    Ab = Ca * bbeta;  % used for computing relative residual norms
    u = b / bbeta;  
    U(:,1) = u;
    
    rb = mvpt(A,u);   % A' * u;  
    
    if tol == 0 && strcmp(type, 'posi')
        r = M \ rb;
    elseif tol == 0 && strcmp(type, 'semi')
        r = Mp * rb;
    else
        % r = pcg(M, rb, tol, 2*n);
        % r = pcg(M, rb, tol, 2*n, M1);
        % r = lsqr([A;L], [u;zeros(p,1)], tol, 2*n);
        r = lsqr(@(z,tflag)afun(z,A,L,tflag),[u;zeros(p,1)],tol,2*n);
    end

    if (isa(A, 'function_handle') || isa(L, 'function_handle'))
        alpha = sqrt(r'*M(r));
        z  = r / alpha;     Z(:,1)  = z;
        Zb(:,1) = M(z);
    else
        alpha = sqrt(r'*M*r);
        z  = r / alpha;     Z(:,1)  = z;
        Zb(:,1) = M * z;
    end
    B(1,1)  = alpha;

    % Prepare for update procedure
    % w = z;
    % phi_bar = bbeta;
    % rho_bar = alpha;
    % x = zeros(n, 1);

    % k step iteration of gGKB and gLSQR
    for j = 1:k 
        fprintf('[gLSQR iterating...], step=%d--------\n', j);
        % compute u in 2-inner product
        s = mvp(A,z) - alpha * u;
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
            fprintf('[Breakdown...], beta=%f, gLSQR breakdown at %d--\n', [beta,j]);
            % U  = U(:,1:j);
            % Z  = Z(:,1:j);
            % Zb = Zb(:,1:j);
            % B  = B(1:j,1:j);
            X  = X(:,1:j);
            res = res(:,1:j);
            bnd = bnd(:,1:j);
            break;
        end

        u = s / beta;
        U(:,j+1) = u;
        B(j+1,j) = beta;

        % compute z in M-inner product
        rb = mvpt(A,u);    % A' * u;
        if tol == 0 && strcmp(type, 'posi')
            r = M \ rb;
        elseif tol == 0 && strcmp(type, 'semi')
            r = Mp * rb;
        else
            % r = pcg(M, rb, tol, 2*n);
            % r = pcg(M, rb, tol, 2*n, M1);
            % r = lsqr([A;L], [u;zeros(p,1)], tol, 2*n);
            r = lsqr(@(z,tflag)afun(z,A,L,tflag),[u;zeros(p,1)],tol,2*n);
        end

        r = r - beta * z;

        if reorth == 1
            for i = 1:j
                r = r - Z(:,i)*(Zb(:,i)'*r);
            end
        elseif reorth == 2
            for i = 1:j
                r = r - Z(:,i)*(Zb(:,i)'*r);
            end
            for i = 1:j
                r = r - Z(:,i)*(Zb(:,i)'*r);
            end
        end
        
        if (isa(A, 'function_handle') || isa(L, 'function_handle'))
            alpha = sqrt(r'*M(r));
        else
            alpha = sqrt(r'*M*r);
        end

        if alpha < 1e-16
            fprintf('[Breakdown...], alpha=%f, gLSQR breakdown at %d--\n', [alpha,j]);
            % U  = U(:,1:j+1);
            % Z  = Z(:,1:j);
            % Zb = Zb(:,1:j);
            % B  = B(1:j+1,1:j);
            X  = X(:,1:j);
            res = res(:,1:j);
            bnd = bnd(:,1:j);
            break;
        end

        z = r / alpha;
        Z(:,j+1) = z;
        if (isa(A, 'function_handle') || isa(L, 'function_handle'))
            Zb(:,j+1)  = M(z);
        else
            Zb(:,j+1)  = M * z;
        end
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
        % rr = A' * (A*x-b);
        % if tol == 0 && strcmp(type, 'posi')
        %     ss = M \ rr;
        % elseif tol == 0 && strcmp(type, 'semi')
        %     ss = Mp * rr;
        %     % ss = M \ rr;
            
        % else
        %     ss = pinv(full(M)) * rr;
        %     % ss = pcg(M, rr, 1e-10, 2*n);
        %     % ss = pcg(M, rr, 1e-4, 2*n, R, R');
        %     % r = lsqr(M, rr,1e-10, 2*n);
        % end

        % rn = sqrt(rr'*ss);
        % res(j) = rn / Ab; 
        % bnd(j) = alpha*beta*abs(yj(j)) / Ab;
        
    end
    
end
    

%--------------------------------------------
function y = afun(z, A, B, transp_flag)
    if strcmp(transp_flag,'transp')   % y = (A(I_n-BB^T))' * z;
        [m,~] = sizemm(A);
        [p,~] = sizemm(B);
        s = mvpt(A, z(1:m));
        t = mvpt(B, z(m+1:m+p));
        y = s + t;
    elseif strcmp(transp_flag,'notransp') % y = (A(I_n-BB^T)) * z;
        s = mvp(A, z);
        t= mvp(B, z);
        y = [s; t];
    end
end