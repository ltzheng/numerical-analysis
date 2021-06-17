clear, clc
format long
h = 10^(0);
val = 1.2;
max_iter = 20;
iter = linspace(1, max_iter, max_iter);
f = @(x) sin(x);
f_prime = @(x) cos(x);
for j = 1:max_iter
    Richardson(j) = N(j, h);
    err(j) = abs(Richardson(j)-f_prime(val));
end
semilogy(iter, err)
xlabel('iter')
ylabel('error')
[min_err, ind] = min(err);
grad = Richardson(ind);
fprintf(['Gradient\t\terror\t\tInitial h\titer.\n']);
fprintf('%.20f\t%d\t%d\t%d\t\n', grad, min_err, h, ind);

function N_j = N(j, h)
    val = 1.2;
    f = @(x) sin(x);
    if j == 1
        N_j = (f(val+h)-f(val))/h;
    else
        item1 = N(j-1,h/2);
        item2 = N(j-1,h);
        N_j = item1+(item1-item2)/(2^(j-1)-1);
    end
end