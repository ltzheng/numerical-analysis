clear, clc
format long
n = 2^6;  % num of points
% define function
f = @(x) sin(2*pi.*x).*exp(1).^(cos(2*pi.*x));
all_x = linspace(0,1,n+1);
all_x = all_x(1:end-1);
y = f(all_x);
% error of lagrange interpolation
test_points = linspace(0,1,1000);
error = zeros(1, 1000);
for i = 1:1000
    u = test_points(i);
    l = zeros(1, n);
    for k = 1:n
        l(k) = l(k)+(-1)^(k-1)/n*sin(n*pi*u)...
            *cot(pi*(u-all_x(k)));
    end
    L = l*y';
    error(i) = abs(L - f(u));
end
% draw semilogy
semilogy(test_points, error)
xlabel('x')
ylabel('Error')
