clear, clc
format long
n = 2^4;  % num of points
% define function
f = @(x) sin(4*(x.^2)) + sin(4*x).^2;
% boundary conditions
m_0 = 0;
m_n = 0;
% cubic spline
all_x = linspace(-1,1,n+1);
y = f(all_x);
h = all_x(2:n+1)-all_x(1:n);
delta = (y(2:n+1)-y(1:n))./h;
d = 6./(h(1:n-1) + h(2:n)).*(delta(2:n)-delta(1:n-1));
lambda = h(2:n)./(h(1:n-1) + h(2:n));
mu = 1.-lambda;
T = diag([mu 1], -1)+diag(2*ones(n+1, 1))...
    +diag([1 lambda], 1);
d_0 = 6/h(1)*(delta(1)-m_0);
d_n = 6/h(n)*(m_n-delta(n));
M = T\[d_0 d d_n]';
% compute error
test_points = linspace(-1,1,2000);
S = [];
for i = 1:2000
    point = test_points(i);
    idx = 0;
    for k = 1:n  % find x_i and x_{i+1}
        if point < all_x(k+1) && point >= all_x(k)
            idx = k;
        end
    end
    if point == 1  % right range value
        idx = n;
    end
    item1 = ((all_x(idx+1)-point)^3*M(idx)...
        +(point-all_x(idx))^3*M(idx+1))/(6*h(idx));
    item2 = ((all_x(idx+1)-point)*y(idx)...
        +(point-all_x(idx))*y(idx+1))/h(idx);
    item3 = h(idx)/6*((all_x(idx+1)-point)*M(idx)...
        +(point-all_x(idx))*M(idx+1));
    S(i) = item1+item2-item3;
end
error = abs(S - f(test_points));
% draw semilogy
semilogy(test_points, error)
xlabel('x')
ylabel('Error')
