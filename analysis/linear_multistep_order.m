clear, clc
format long
f = @(x, y) x*exp(-4*x)-4*y;
g = @(x) x*exp(-4*x);
acc = @(x) 1/2*x^2*exp(-4*x);
left = 0;
right = 2;
val = 2;
steps = [];
err = [];
for m = 10:10:500
    x = linspace(left, right, m);
    h = (right-left)/m;
    steps = [steps, h];
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
    err = [err, abs(acc(val)-y(m)-(acc(val-h)-y(m-1)))];
end
loglog(steps, err)
xlabel('steps')
ylabel('error')