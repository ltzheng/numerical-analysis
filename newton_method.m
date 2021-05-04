clear, clc

syms x
f(x) = x^3 - 3*x^2 + 2;
f_prime(x) = diff(f(x));
g(x) = x - f(x) / f_prime(x);

% initial settings
x_0 = [-2.5, 0.5, 2.5];
tolerance = 1e-15;

N = 300;

x_l = newton(x_0(1), g, tolerance, N);
x_m = newton(x_0(2), g, tolerance, N);
x_r = newton(x_0(3), g, tolerance, N);

order_l = newton_with_order(x_0(1), x_l, g, tolerance, N);
order_m = newton_with_order(x_0(2), x_m, g, tolerance, N);
order_r = newton_with_order(x_0(3), x_r, g, tolerance, N);

function root = newton(x_0, g, tolerance, N)
    for iter = 1:N
        x_1 = vpa(g(x_0));
        if abs(x_1 - x_0) < tolerance
            root = vpa(x_1);
            break
        end
        x_0 = x_1;
    end
end

function orders = newton_with_order(x_0, x_exact, g, tolerance, N)
    fprintf('initial: x = %f\n', x_0)
    fprintf('iter.\tx\t\torder\n');
    dashString = repmat('-', 1, 30);
    fprintf('%s\n', dashString);

    x_1 = vpa(g(x_0));
    e_kp1 = abs(x_1 - x_exact);
    e_k = abs(x_1 - x_exact);
    e_km1 = abs(x_1 - x_exact);

    for iter = 1:N
        if abs(x_1 - x_0) < tolerance
            root = vpa(x_1);
            fprintf('%s\n', dashString);
            fprintf('True solution: x = %f\n\n', root)
            break
        end
        if iter >= 3
            orders(iter, 1) = converge_order(e_kp1, e_k, e_km1);
            fprintf('%d\t%f\t%f\n', iter, x_0, orders(iter, 1));
        else
            fprintf('%d\t%f\n', iter, x_0);
        end
        x_0 = x_1;
        e_km1 = e_k;
        e_k = e_kp1;
        x_1 = vpa(g(x_0));
        e_kp1 = abs(x_1 - x_exact);
    end
end

function order = converge_order(e_kp1, e_k, e_km1)
    order = log(abs(e_kp1 / e_k)) / log(abs(e_k / e_km1));
end
