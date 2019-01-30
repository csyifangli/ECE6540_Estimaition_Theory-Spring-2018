A=[1 1;1 1];
b=[0.5 0.5;0.6 0.6];
gamma=100;
    cvx_begin
        variable x(2,2);
        x_1=x(1,1)+x(1,2);
        x_2=x(2,1)+x(2,2);
        y=[x_1^2,x_2^2]';
        minimize( norm(A*x-b)+gamma*norm(y,1));
    cvx_end
x