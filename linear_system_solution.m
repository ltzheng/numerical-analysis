clear, clc

SIZE = 10;

A = 2 * eye(SIZE);
for i = 1:SIZE
    if i ~= 1 && i ~= SIZE
        A(i, i - 1) = -1;
        A(i, i + 1) = -1;
    end
    if i == 1
        A(i, i + 1) = -1;
    end
    if i == SIZE
        A(i, i - 1) = -1;
    end
end


b = [2 -2 2 -1 0 0 1 -2 2 -2]';

N = 500;

global exact
exact = [1 0 1 0 0 0 0 -1 0 -1]';

[jacobi_solution, jacobi_error] = jacobi(A, b, N);
[gs_solution, gs_error] = gauss_seidel(A, b, N);

jacobi_solution
gs_solution

x_range = 1:N;
semilogy(x_range, jacobi_error, '-', x_range, gs_error, '--')
xlabel('Iteration')
ylabel('Error')
hold on

omegas = [0.2, 0.8, 1.2, 1.4, 1.6, 1.7];
for i = 1:size(omegas, 2)
    [sor_solution, sor_error] = SOR(A, b, N, omegas(i));
    semilogy(x_range, sor_error)
    legend_str{i} = ['\omega=' num2str(omegas(i))];
end

legend(['Jacobi', 'Gauss-Seidel', legend_str])

grid on

function [solution, error] = jacobi(A, b, N)
    global exact
    error = zeros(N, 1);

    SIZE = size(b, 1);
    
    D = diag(diag(A));
    inversed_D = zeros(SIZE, SIZE);
    for i = 1:SIZE
        inversed_D(i, i) = 1 / D(i, i);
    end
    R = eye(SIZE) - inversed_D * A;
    g = inversed_D * b;

    solution = zeros(SIZE, 1);
    
    % iteration
    for iter = 1:N
        solution = R * solution + g;

        % compute error
        e = norm(solution - exact, inf);
        error(iter, 1) = e;
    end
end

function [solution, error] = gauss_seidel(A, b, N)
    global exact
    error = zeros(N, 1);

    SIZE = size(b, 1);
    
    D = diag(diag(A));
    L = tril(A) - D;
    U = triu(A) - D;
    inversed_DplusL = (D + L) \ eye(SIZE);
    S = -inversed_DplusL * U;
    f = inversed_DplusL * b;

    solution = zeros(SIZE, 1);
    
    % iteration
    for iter = 1:N
        solution = S * solution + f;
        
        % compute error
        e = norm(solution - exact, inf);
        error(iter, 1) = e;
    end
end

function [solution, error] = SOR(A, b, N, omega)
    global exact
    error = zeros(N, 1);

    SIZE = size(b, 1);
    
    D = diag(diag(A));
    L = tril(A) - D;
    U = triu(A) - D;

    inversed_D = D \ eye(SIZE);
    inversed_DplusL = (eye(SIZE) + omega * inversed_D * L) \ eye(SIZE);
    S = inversed_DplusL * ((1 - omega) * eye(SIZE) - omega * inversed_D * U);
    f = omega * inversed_DplusL * inversed_D * b;

    solution = zeros(SIZE, 1);
    
    % iteration
    for iter = 1:N
        solution = S * solution + f;
        
        % compute error
        e = norm(solution - exact, inf);
        error(iter, 1) = e;
    end
end