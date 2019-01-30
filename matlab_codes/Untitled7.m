gamma = logspace( -2,-0.1, 10 );
l2norm = zeros(size(gamma));
l1norm = zeros(size(gamma));
fprintf( 1, '   gamma       norm(x,1)    norm(A*x-b)\n' );
fprintf( 1, '---------------------------------------\n' );
A = rand(4,4);
b=[exp(1+j) 1 1 1]';
for k = 1:length(gamma),
    fprintf( 1, '%8.4e', gamma(k) );
    cvx_begin
        variable x(4);
        minimize( norm(A*x-b)+gamma(k)*norm(x,1) );
    cvx_end
    l1norm(k) = norm(x,1);
    l2norm(k) = norm(A*x-b);
    fprintf( 1, '   %8.4e   %8.4e\n', l1norm(k), l2norm(k) );
end
plot( l1norm, l2norm );
xlabel( 'norm(x,1)' );
ylabel( 'norm(A*x-b)' );
grid on
x