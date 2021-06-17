clear, clc
format long
max_range_num = 20;
m_s = linspace(1, max_range_num, max_range_num);
f = @(x) exp(cos(x));
% calculated by Wolfram alpha
true_integral = 7.954926521012845274513219665329394328161342771816638573400595955383360608164694666995137357228568774;
left = -pi;
right = pi;
for m = 1:max_range_num
    h = (right-left)/m;
    item = 1/2*f(left)+1/2*f(right);
    for i = 1:m-1
        item = item+f(left+i*h);
    end
    trapezoid(m) = h*item;
    err(m) = abs(trapezoid(m)-true_integral);
end
semilogy(m_s, err)
xlabel('m')
ylabel('error')
[min_err, ind] = min(err);
result = trapezoid(ind);
fprintf(['result\t\terror\tm\n']);
fprintf('%.20f\t%2d\t%d\n', result, min_err, ind);
