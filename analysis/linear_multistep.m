clear, clc
format long
f = @(x, y) x*exp(-4*x)-4*y;
g = @(x) x*exp(-4*x);
m = 500;
left = 0;
right = 2;
h = (right-left)/m;
x = linspace(left, right, m);
% 2-order Runge-Kutta
y(1) = 0;
k1 = f(x(1), y(1));
k2 = f(x(1)+h/2, y(1)+h/2*k1);
y(2) = y(1)+h*k2;
% linear multistep
for n = 2:m-1
    y(n+1) = (1/(1+4/3*h))* ...
        (1/3*h*g(x(n+1))+y(n-1) ...
        +4/3*h*f(x(n),y(n))+ ...
        1/3*h*f(x(n-1),y(n-1)));
end
plot(x, y)
xlabel('x')
ylabel('y')