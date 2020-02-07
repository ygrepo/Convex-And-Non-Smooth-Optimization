function [x, y] = gen_random_vectors(n, p)
    if (isinf(p))
        pnorm = 1;
    else
        pnorm = p;
    end
    r = randn(2,n); % Use a large n
    norm_r = (sum(abs(r).^pnorm,1)).^(1/pnorm);
    r = r./norm_r;
    %r = bsxfun(@rdivide, r, norm_r);
    x = r(1,:);
    y = r(2,:);
