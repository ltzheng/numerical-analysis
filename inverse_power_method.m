clear, clc

rng(2);
A = rand(100, 100);
target = 0.8 - 0.6i;

tolerance = 1e-13;
N = 1000;
inv_shift_power_m(A, N, target, tolerance)

function inv_shift_power_m(A, N, p, tolerance)
    A = A - p * eye(size(A, 1));  % shift
    q_old = ones(size(A, 1), 1);
    q_old_bar = q_old / max(q_old);
    q_old_bar_gap = q_old_bar;
    q_new_bar_gap = q_old_bar;
    for iter = 1:N
        q_new = A \ q_old_bar;  % inverse power method
        mu = max(q_new);  % eigenvalue
        q_new_bar = q_new / mu;  % eigenvector
        % only one dominant eigenvalue
        if norm(q_old_bar - q_new_bar, inf) < tolerance
            fprintf('positive dominant eigenvalue:');
            lambda = p + 1 / mu
            q_new_bar
            break
        end
        if norm(q_old_bar + q_new_bar, inf) < tolerance
            fprintf('negative dominant eigenvalue');
            lambda = -lambda
            q_new_bar
            break
        end
        % two opposite dominant eigenvalues
        if iter > 3
            if (norm(q_old_bar_gap - q_old_bar, inf) < tolerance) ...,
                && (norm(q_new_bar_gap - q_new_bar, inf) < tolerance)
                fprintf('two opposite dominant eigenvalues');
                % consider A^2, yield \lambda^2
                q_old = q_new;
                q_new = A \ q_new;
                mu = sqrt(q_new(1) / q_old_bar(1));
                lambda1 = p + 1 / mu
                lambda2 = -lambda1
                q_new_bar1 = q_new + lambda1 * q_old;
                q_new_bar2 = q_new + lambda2 * q_old;
                q_new_bar1 = q_new_bar1 / max(q_new_bar1)
                q_new_bar2 = q_new_bar2 / max(q_new_bar2)
                break
            end
        end
        q_old_bar_gap = q_new_bar_gap;  % update recent 4 q-vectors
        q_new_bar_gap = q_old_bar;
        q_old_bar = q_new_bar;
    end
end
