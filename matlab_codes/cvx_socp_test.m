%sample cvx SOCP
f=[1 2]';
A_1=[1 -1;1 1];
A_2=[0 1;1 1];
b_1=[1 0]';
b_2=[0.5 1]';
c_1=[1 0]';
c_2=[0.5 0.5]';

cvx_begin
variable x(2)
minimize (f'*x);
subject to
norm(A_1*x+b_1)<= 10;
%norm(A_2*x+b_2)<= c_2'*x;
x(1)<=x(2)-6;
cvx_end
x
