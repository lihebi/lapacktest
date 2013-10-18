function final(n)
    A=randn(n,n); b=randn(n,1); Abk=A; pvt = 1:n;
    A1=A;pvt1=pvt;
    [A, pvt] = mydgetrf(A1, pvt1)
    x = mydtrsm(A, b, pvt);
    %can also use:
    % y = mydtrsm_f(A, b, pvt)
    % x = mydtrsm_b(A, y)
    xx= Abk\b;
    norm(x'-xx)
end