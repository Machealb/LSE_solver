function [x, x1, x2, b, d] = gen_prob2(A, C, type)
    % construct test problem for Linear least squares problem 
    % with linear equaily constraints:
    %    min ||Ax-b||_2   s.t.   ||Cx-d||_2=min.
    % 
    % Inputs:
    %   A, C: the matrices of LSE
    %   type: control the type of true (discretized) solution by setting w
    %       '1': w = f(t) = t, t \in [0,1]       
    %       '2': w = f(t) = t-t^2,  t \in [0,1]
    % 
    % Outputs:
    %   x: true solution
    %   b, d: right-hand terms
    %
    % Haibo Li, School of Mathematics and Statistics, The University of Melbourne
    % 31, Dec, 2024.

    [~, n] = size(A);
    if size(C,2) ~= n
        error('The dimensions are not consistent')
    end

    if strcmp(type, '1')
        tt = linspace(0,1,n);
        f = @(t) t;
        x1 = f(tt);
        x1 = x1(:);
    elseif strcmp(type, '2')
        tt = linspace(-1,1,n);
        f = @(t) t.^3 + t.^2;
        x1 = f(tt);
        x1 = x1(:);
    elseif strcmp(type, '3')
        tt = linspace(-pi,pi,n);
        f = @(t) sin(5*t)+2*cos(t);
        x1 = f(tt);
        x1 = x1(:);
    end

    % construct x1
    G = A'*A + C'*C;
    V = null(full(C));
    B = V' * G * V;
    x1 = x1 - V * (B \ (V'*(G*x1)));    % x'*G*N(C)=0

    d1 = C * x1;
    % res1 = randn(m,1);
    % res1 = res1 / norm(res1);
    % res  = res1 - A * (A \ res1);
    % b = b1 + res/norm(res);
    d = d1;

    % construct x2
    M = (A*V)';
    w2 = 10*M(:,1);
    x2 = V * w2;
    b  = A * x2;

    x = x1 + x2;

end