clear, clc
format long
n = 2^4;  % num of points
% define function
syms x;
f(x) = sin(4*(x.^2)) + sin(4*x).^2;
f_prime(x) = diff(f(x));
% boundary conditions
m_0 = vpa(f_prime(-1));
m_n = vpa(f_prime(1));
% cubic spline
all_x = linspace(-1,1,n);
y = vpa(f(all_x));
h = all_x(2:n)-all_x(1:n-1);
delta = (y(2:n)-y(1:n-1))./h;
d = 6./(h(1:end-1) + h(2:end)).*(delta(2:end)-delta(1:end-1));
lambda = h(1:end-1)./(h(1:end-1) + h(2:end));
mu = 1.-lambda;
T = diag([mu 1], -1)+diag(2*ones(n, 1))+diag([1 lambda], 1);
d_0 = 6./h(1)*(delta(1)-m_0);
d_n = 6./h(n-1)*(m_n-delta(n-1));
M = T\[d_0 d d_n]'
% compute error
error = zeros(2000, 1);
test_points = linspace(-1,1,2000);
for i = 1:2000
    point = test_points(i);
    idx = min(find(all_x<=point,1,'last'), n-1);
    x_i = all_x(idx);
    x_ip1 = all_x(idx+1);
    M_i = M(idx);
    M_ip1 = M(idx+1);
    y_i = y(idx);
    y_ip1 = y(idx+1);
    h_i = h(idx);
    item1 = ((x_ip1-point)^3*M_i+(point-x_i)^3*M_ip1)/(6*h_i);
    item2 = ((x_ip1-point)*y_i+(point-x_i)*y_ip1)/h_i;
    item3 = h_i/6*((x_ip1-point)*M_i+(point-x_i)*M_ip1);
    S = item1+item2-item3;
    error(i, 1) = abs(S - vpa(f(point)));
end
% draw semilogy
semilogy(test_points, error)
xlabel('x')
ylabel('Error')
% hold on
% % draw loglog