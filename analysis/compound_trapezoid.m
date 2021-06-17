clear, clc
format long
max_range_num = 30;
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

% clear ,clc;
% % define function and pre - define real value
% f = @(x)exp(cos(x));
% real_val = 7.954926521012845274513219665329394328;
% for k = 1 : 18
% x(k) = k;
% error (k) = abs( real_val - ...
% my_integrate (f, x(k), -pi , pi ));
% end
% semilogy (x, error );
% % calculate integrate value
% function ans = my_integrate (f, m, start , end_)
% ans = 0;
% h = (end_ - start ) / m;
% ans = ans + h * (f( start ) + f( end_ )) / 2;
% for k = 1 : m - 1
% ans = ans + h * f( start + k * h);
% end
% end