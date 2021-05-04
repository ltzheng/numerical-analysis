clear, clc

% A = [-148 -105 -83 -67; 488 343 269 216; -382 -268 -210 -170; 50 38 32 39];
% A = -A;
A = [222 580 584 786; -82 -221 -208 -288; 37 98 101 132; -30 -82 -88 -109];

tolerance = 1e-12;
N = 500;
[lambda, q_new_bar] = power_m(A, N, tolerance)

function [lambda, q_new_bar] = power_m(A, N, tolerance)
    q_old = ones(size(A, 1), 1);
    q_old_bar = q_old / norm(q_old, inf);
    for iter = 1:N
        q_new = A * q_old_bar;
        lambda = norm(q_new, inf);  % eigen value
        q_new_bar = q_new / lambda;  % eigen vector

        % only one dominant eigenvalue
        pos = norm(q_old_bar - q_new_bar, inf) < tolerance;
        neg = norm(q_old_bar + q_new_bar, inf) < tolerance;
        if pos || neg
            if pos
                fprintf('positive eigenvalue');
                break
            end
            if neg
                fprintf('negative eigenvalue');
                lambda = -lambda;
                break
            end
        end

        % two opposite eigenvalues
        if iter > 3
            sep_converge = (norm(q_old_bar_gap - q_old_bar, inf) < tolerance) && (norm(q_new_bar_gap - q_new_bar, inf) < tolerance);
            if sep_converge
                q_old_bar = q_new_bar;
                q_new_bar = A * q_new_bar;  % consider A^2
                squared_lambda = norm(q_new_bar, inf) / norm(q_old_bar, inf);
                lambda1 = sqrt(squared_lambda);
                lambda2 = -lambda1;
                fprintf('two opposite eigenvalue');
                lambda = [lambda1, lambda2];
                q_new_bar1 = q_new_bar + lambda_1 * q_old_bar;
                q_new_bar2 = q_new_bar + lambda_2 * q_old_bar;
                q_new_bar = [q_new_bar1, q_new_bar2];
            end
        end
        q_new_bar_gap = q_old_bar;
        q_old_bar_gap = q_new_bar_gap;
        q_old_bar = q_new_bar;
    end
end