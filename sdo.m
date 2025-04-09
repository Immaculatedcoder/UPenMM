function D2 = sdo(n)
    h = 1/n;
    e = ones(n,1);
    D2 = spdiags([e -2*e e], -1:1, n, n);
    D2(1, n) = 1;
    D2(n, 1) = 1;

    % Scale by h^2
    D2 = D2 / h^2;
end