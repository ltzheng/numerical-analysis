clear, clc
format long
n_s = [2^4, 2^5, 2^6, 2^7, 2^8, 2^9, 2^10];% num of points
errors = zeros(1, size(n_s, 2));
for i = 1:size(n_s, 2)
    errors(i) = max(cubic(n_s(i)));
end
% draw loglog
loglog(n_s, errors)
xlabel('n')
ylabel('Max error')

function error = cubic(n)
    % define function
    f = @(x) sin(4*(x.^2)) + sin(4*x).^2;
    % cubic spline
    all_x = linspace(-1,1,n+1);
    y = f(all_x);
    h = all_x(2:n+1)-all_x(1:n);
    delta = (y(2:n+1)-y(1:n))./h;
    d = 6./(h(1:n-1) + h(2:n)).*(delta(2:n)-delta(1:n-1));
    lambda = h(2:n)./(h(1:n-1) + h(2:n));
    mu = 1.-lambda;
    T = diag([mu 1], -1)+diag(2*ones(n+1, 1))...
    +diag([0 lambda], 1);
    T(1, n+1) = -1;
    T(n+1, 1) = h(1)/6;
    T(n+1, 2) = h(1)/3;
    T(n+1, n) = h(n)/3;
    T(n+1, n+1) = h(n)/6;
    % boundary conditions
    d_0 = 0;
    d_n = delta(n)-delta(1);
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
end
