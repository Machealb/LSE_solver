function [A, C, x, b, d, x1, x2] = gen_prob3(n)
    % construct test problem for Linear least squares problem 
    % with linear equaily constraints:
    %    min ||Ax-b||_2   s.t.   ||Cx-d||_2=min.
    % 
    % Inputs:
    %   A_type: for generating A
    %   C_type: for generating C=L
    %   x_type: control the type of true (discretized) solution by setting w
    %       '1': w = f(t) = t, t \in [0,1]       
    %       '2': w = f(t) = t-t^2,  t \in [0,1]
    % 
    % Outputs:
    %   A, C: matrices
    %   x: true solution
    %   b, d: right-hand terms
    %
    % Haibo Li, School of Mathematics and Statistics, The University of Melbourne
    % 30, May, 2025.

    % if strcmp(A_type, '1')
    %     A = get_l(n,1);
    %     % V = ones(n,1);
    % elseif strcmp(A_type, '2')
    %     A = get_l(n,2);
    %     % v1 = ones(n,1);
    %     % v2 = 1:n;
    %     % v2 = v2(:);
    %     % V = [v1,v2];
    % elseif strcmp(A_type, '3')
    %     A = mmread('lp_bnl2.mtx');
    %     n = size(A,2);
    % end

    % if strcmp(C_type, '1')
    %     r0 = 160;
    %     r = n - r0;
    %     c = zeros(n,1);
    %     c(1:r) = linspace(10,0.01,r);
    %     C = spdiags(c(:),0,n,n); 
    %     E = speye(n);
    %     V = E(:,r+1:n);
    % else
    %     error('no setting yet')
    % end


    r1 = 200;
    r2 = 300;
    r3 = n - r1 - r2;
    c = zeros(n,1);
    c(1:r1) = ones(r1,1);
    c(r1+1:r1+r2) = linspace(0.99,0.01,r2);
    s = zeros(n,1);
    cr = c(r1+1:r1+r2);
    s(r1+1:r1+r2) = sqrt(1 - cr.^2);
    s(r1+r2+1:n) = ones(r3,1);
    C1 = spdiags(c(:),0,n,n);
    S1 = spdiags(s(:),0,n,n);
    d1 = linspace(1,100,n);
    D  = spdiags(d1(:),0,n,n);
    invD = spdiags(1.0./d1(:),0,n,n);
    A = C1 * D;
    C = S1 * D;
    

    % construct x1
    E = speye(n);
    V = invD * E(:,1:r1);
    G = A'*A + C'*C;
    z = [zeros(r1,1);ones(r2,1);zeros(r3,1)];
    x1 = G * z;

    d1 = C * x1;
    % res1 = randn(m,1);
    % res1 = res1 / norm(res1);
    % res  = res1 - A * (A \ res1);
    % b = b1 + res/norm(res);
    d = d1;

    % construct x2
    % M = (A*V)';
    % w2 = M(:,1);
    ww = linspace(100,1,r1);
    w2 = ww(:);
    x2 = V * w2;
    b  = A * x2;

    x = x1 + x2;

end
