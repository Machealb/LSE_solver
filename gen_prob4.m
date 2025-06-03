function [Afun, Cfun, x, b, d, x1, x2] = gen_prob4(n, t, x_type)
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

    Afun = diff2_operator(n);
    [Cfun, V] = structured_operator(n, t);

    if strcmp(x_type, '1')
        tt = linspace(0,1,n);
        f = @(t) t;
        x1 = f(tt);
        x1 = x1(:);
    elseif strcmp(x_type, '2')
        tt = linspace(-1,1,n);
        f = @(t) t.^2;
        x1 = f(tt);
        x1 = x1(:);
    elseif strcmp(x_type, '3')
        tt = linspace(-pi,pi,n);
        f = @(t) sin(2*t)+3*cos(t);
        x1 = f(tt);
        x1 = x1(:);
    elseif strcmp(x_type, '4')
        x1 = ones(n,1);
    end

    % construct x1
    Gfun = @(x) Afun(Afun(x,'notransp'),'transp') + Cfun(Cfun(x,'notransp'),'transp');
    B = V' * Gfun(V);
    x1 = x1 - V * (B \ (V'*Gfun(x1)));    % x'*G*N(C)=0
    d1 = Cfun(x1,'notransp');
    % [p,~] = sizemm(Cfun);
    % res1 = ones(p,1);
    % res1 = res1 / norm(res1);
    % d = d1 + res1;
    d = d1;

    % construct x2
    M1 = Afun(V,'notransp');
    M = M1';
    ll = linspace(100,10,400);
    w2 = M(:,1:400)*ll';
    x2 = V * w2;
    b  = Afun(x2,'notransp');

    x = x1 + x2;

end
