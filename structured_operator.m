function [Cfun, N] = structured_operator(n, d)
% MAKE_STRUCTURED_OPERATOR_C constructs a matrix-free operator C = M * P'
% where P is a built-in orthogonal matrix from gallery('orthog', n).
%
% Inputs:
%   n - dimension of the domain (C: R^n-->R^r)
%   d - desired null space dimension (must satisfy d <= 1000)
%
% Outputs:
%   Cfun - function handle for matrix-vector products with C and C^T
%   N    - explicit basis for the null space of C

    assert(d <= 1000, 'Null space dimension d must be <= 1000.');
    r = n - d;

    % --- Step 1: Generate orthogonal matrix using built-in method
    P = gallery('orthog', n, 1);  % Fast and deterministic orthogonal matrix
    % ll = linspace(2,1,r);
    ll = ones(r,1);
    ll = ll(:);

    % --- Step 2: Define matrix-free operator
    Cfun = @(x, mode) apply_C(x, P, ll, r, d, mode);

    % --- Step 3: Construct null space basis N ∈ R^{n × d}
    N = P(:, r+1:end);
end

function y = apply_C(x, P, ll, r, d, mode)
% APPLY_C evaluates C*x or C'*x where C = M * P', M = [I_r  0]

    n = size(P, 1);  % just in case
    if isempty(x) && strcmp(mode, 'size')
        % Return size of the operator: [rows, cols]
        y = [r, n];
        return;
    end

    switch mode
        case 'notransp'
            % Forward: C * x = (P' * x)(1:r)
            y = P' * x;
            y = diag(ll) * y(1:r,:);

        case 'transp'
            % Transpose: C' * y = P * [y; 0]
            l = size(x,2);
            y = P * [diag(ll)*x; zeros(d,l)];

        otherwise
            error('Invalid mode. Use ''notransp'' or ''transp''.');
    end
end
