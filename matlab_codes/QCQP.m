%========QCQP opt. solving======================
% by Xiang Zhang
m = 16; n = 8;
A = randn(m,n);
b = randn(m,1);

cvx_begin
    variable x(n)
    minimize( norm(A*x-b) )
cvx_end
x
x_ls = A \ b