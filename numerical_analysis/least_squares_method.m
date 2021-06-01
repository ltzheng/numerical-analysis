clear, clc
format long
syms a b
x_i = [2.1,2.5,2.8,3.2];
y_i = [0.6087,0.6849,0.7368,0.8111];
Q = sum(((a+b*x_i)./x_i-1./y_i).^2);
dQa = diff(Q, a);
dQb = diff(Q, b);
dQa_coef = coeffs(dQa);
dQb_coef = coeffs(dQb);
coef=[dQa_coef(3),dQa_coef(2);dQb_coef(3),dQb_coef(2)]
y=-[dQa_coef(1);dQb_coef(1)]
ab=coef\y;
a=ab(1)
b=ab(2)
f = @(x) x./(a+b*x);
pred = f(x_i);
error = double(norm(pred-y_i, 2))
