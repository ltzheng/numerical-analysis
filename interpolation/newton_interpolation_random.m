clear, clc
format long
rng(22);
n_s = [2^2, 2^3, 2^4, 2^5, 2^6, 2^7];% num of points
errors = zeros(1, size(n_s, 2));
for i = 1:size(n_s, 2)
    errors(i) = max(newton(n_s(i)));
end
% draw semilogy
semilogy(n_s, errors)
xlabel('n')
ylabel('Max error')

function error = newton(n)
    % define function
    f = @(x) 1./(1+25*x.^2);
    p = randperm(n+1);
    for i = 1:n+1  % chebyshev
        all_x(i) = cos((p(i)-1)*pi/n);
    end
    % table of difference quotient
    g = f(all_x);
    for k = 1:n  % compute k-order diff
        for j = n+1:-1:k+1
            g(j) = (g(j)-g(j-1))/(all_x(j)-all_x(j-k));
        end
    end
    % error of newton interpolation
    test_points = linspace(-1,1,2000);
    error = zeros(1, 2000);
    for i = 1:2000
        u = test_points(i);
        t = 1;
        val = g(1);
        for k = 2:n+1
            t = t*(u-all_x(k-1));
            val = val+t*g(k);
        end
        error(i) = abs(val - f(u));
    end
end
