clear, clc

% A = [-148 -105 -83 -67; 488 343 269 216; -382 -268 -210 -170; 50 38 32 29];
% A = -A;
A = [222 580 584 786; -82 -211 -208 -288; 37 98 101 132; -30 -82 -88 -109];

tolerance = 1e-15;
N = 1000;
power_m(A, N, tolerance);

function power_m(A, N, tolerance)
    q_old = ones(size(A, 1), 1);
    q_old_bar = q_old / norm(q_old, inf);
    q_old_bar_gap = q_old_bar;
    q_new_bar_gap = q_old_bar;
    for iter = 1:N
        q_new = A * q_old_bar;
        lambda = norm(q_new, inf);  % eigenvalue
        q_new_bar = q_new / lambda;  % eigenvector

        % only one dominant eigenvalue
        if norm(q_old_bar - q_new_bar, inf) < tolerance
            fprintf('positive dominant eigenvalue:');
            lambda
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
                q_new = A * q_new;
                lambda1 = sqrt(q_new(1) / q_old_bar(1))
                lambda2 = -lambda1
                q_new_bar1 = q_new + lambda1 * q_old;
                q_new_bar2 = q_new + lambda2 * q_old;
                q_new_bar1 = q_new_bar1 / norm(q_new_bar1, inf)
                q_new_bar2 = q_new_bar2 / norm(q_new_bar2, inf)
                break
            end
        end
        q_old_bar_gap = q_new_bar_gap;  % update recent 4 q-vectors
        q_new_bar_gap = q_old_bar;
        q_old_bar = q_new_bar;
    end
end