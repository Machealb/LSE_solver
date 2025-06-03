function Afun = diff2_operator(n)
% Creates second-order finite difference operator A ∈ ℝ^{(n-2) × n}
% Usage: Afun = make_diff2_operator(n);
%        y = Afun(x, 'notransp');   % A * x
%        y = Afun(x, 'transp');     % A' * x
    Afun = @(x, mode) apply_diff2_operator(x, n, mode);
end

function y = apply_diff2_operator(x, n, mode)
    if isempty(x) && strcmp(mode, 'size')
        % Return size of the operator: [rows, cols]
        y = [n-2, n];
        return;
    end

    if strcmp(mode, 'notransp')
        y = x(1:end-2,:) - 5*x(2:end-1,:) + x(3:end,:);
        y = 0.1 * y;
    elseif strcmp(mode, 'transp')
        l = size(x,2);
        y = zeros(n, l);
        y(1,:)   = x(1,:);
        y(2,:)   = -5*x(1,:) + x(2,:);
        y(3:n-2,:) = x(1:end-2,:) - 5*x(2:end-1,:) + x(3:end,:);
        y(n-1,:) = -5*x(end-1,:) + x(end,:);
        y(n,:)   = x(end,:);
        y = 0.1 * y;
    else
        error('Mode must be ''notransp'' or ''transp''.');
    end
end
