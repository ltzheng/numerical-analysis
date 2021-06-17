clear, clc
format long
h_exp = linspace(0,-15,16);
h = 10.^h_exp;
val = 1.2;
f = @(x) sin(x);
f_prime = @(x) cos(x);
difference = (f(val+h)-f(val))./h;
err = abs(difference-f_prime(val));
loglog(h, err)
xlabel('h')
ylabel('error')