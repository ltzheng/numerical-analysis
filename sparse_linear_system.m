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
tolerance = 1e-15;
N = 2000;
omega = 1.6;
TEST_NUM = 100;

global exact
exact = [1 0 1 0 0 0 0 -1 0 -1]';

fprintf('Orgiginal Jacobi')
tic
for i = 1:TEST_NUM
    [jacobi_solution, jacobi_error] = jacobi(A, b, N, tolerance);
end
toc
fprintf('Orgiginal Gauss-Seidel')
tic
for i = 1:TEST_NUM
    [gs_solution, gs_error] = gauss_seidel(A, b, N, tolerance);
end
toc
fprintf('Orgiginal SOR')
tic
for i = 1:TEST_NUM
    [sor_solution, sor_error] = SOR(A, b, N, omega, tolerance);
end
toc

fprintf('Optimized Jacobi')
tic
for i = 1:TEST_NUM
    [jacobi_solution, jacobi_error] = sparse_jacobi(A, b, N, tolerance);
end
toc
fprintf('Optimized Gauss-Seidel')
tic
for i = 1:TEST_NUM
    [gs_solution, gs_error] = sparse_gauss_seidel(A, b, N, tolerance);
end
toc
fprintf('Optimized SOR')
tic
for i = 1:TEST_NUM
    [sor_solution, sor_error] = sparse_SOR(A, b, N, omega, tolerance);
end
toc

function [x2, error] = jacobi(A, b, N, tolerance)
    global exact
    error = zeros(N, 1);

    SIZE = size(b, 1);
    
    x1 = zeros(SIZE, 1);
    x2 = ones(SIZE, 1);
    
    % iteration
    for iter = 1:N
        x1 = x2;
        for i = 1:SIZE
            s = 0;
            for j = 1:SIZE
                s = s + A(i, j) * x1(j);
            end
            x2(i) = (b(i) - s + A(i, i) * x1(i)) / A(i, i);
        end
        % compute error
        e = norm(x2 - exact, inf);
        error(iter, 1) = e;
        if e < tolerance
            break
        end
    end
end

function [x2, error] = gauss_seidel(A, b, N, tolerance)
    global exact
    error = zeros(N, 1);

    SIZE = size(b, 1);
    x1 = zeros(SIZE, 1);
    x2 = ones(SIZE, 1);
    
    % iteration
    for iter = 1:N
        x1 = x2;
        for i = 1:SIZE
            s = 0;
            for j = 1:SIZE
                s = s + A(i, j) * x2(j);
            end
            x2(i) = (b(i) - s + A(i, i) * x2(i)) / A(i, i);
        end
        % compute error
        e = norm(x2 - exact, inf);
        error(iter, 1) = e;
        if e < tolerance
            break
        end
    end
end

function [x2, error] = SOR(A, b, N, omega, tolerance)
    global exact
    error = zeros(N, 1);

    SIZE = size(b, 1);
    x1 = zeros(SIZE, 1);
    x2 = ones(SIZE, 1);
    
    % iteration
    for iter = 1:N
        x1 = x2;
        for i = 1:SIZE
            s = 0;
            for j = 1:SIZE
                if i ~= j
                    s = s + A(i, j) * x2(j);
                end
            end
            x2(i) = (omega * (b(i) - s) + (1 - omega) * A(i, i) * x2(i)) / A(i, i);
        end

        % compute error
        e = norm(x2 - exact, inf);
        error(iter, 1) = e;
        if e < tolerance
            break
        end
    end
end

function [x2, error] = sparse_jacobi(A, b, N, tolerance)
    global exact
    error = zeros(N, 1);

    SIZE = size(b, 1);
    
    x1 = zeros(SIZE, 1);
    x2 = ones(SIZE, 1);
    
    % iteration
    for iter = 1:N
        x1 = x2;
        for i = 1:SIZE
            s = 0;
            if i ~= 1 && i ~= SIZE
                low = i - 1;
                high = i + 1;
            else
                if i == 1
                    low = i;
                    high = i + 1;
                end
                if i == SIZE
                    low = i - 1;
                    high = i;
                end
            end
            for j = low:high
                s = s + A(i, j) * x1(j);
            end
            x2(i) = (b(i) - s + A(i, i) * x1(i)) / A(i, i);
        end
        % compute error
        e = norm(x2 - exact, inf);
        error(iter, 1) = e;
        if e < tolerance
            break
        end
    end
end

function [x2, error] = sparse_gauss_seidel(A, b, N, tolerance)
    global exact
    error = zeros(N, 1);

    SIZE = size(b, 1);
    x1 = zeros(SIZE, 1);
    x2 = ones(SIZE, 1);
    
    % iteration
    for iter = 1:N
        x1 = x2;
        for i = 1:SIZE
            s = 0;
            if i ~= 1 && i ~= SIZE
                low = i - 1;
                high = i + 1;
            else
                if i == 1
                    low = i;
                    high = i + 1;
                end
                if i == SIZE
                    low = i - 1;
                    high = i;
                end
            end
            for j = low:high
                s = s + A(i, j) * x2(j);
            end
            x2(i) = (b(i) - s + A(i, i) * x2(i)) / A(i, i);
        end
        % compute error
        e = norm(x2 - exact, inf);
        error(iter, 1) = e;
        if e < tolerance
            break
        end
    end
end

function [x2, error] = sparse_SOR(A, b, N, omega, tolerance)
    global exact
    error = zeros(N, 1);

    SIZE = size(b, 1);
    x1 = zeros(SIZE, 1);
    x2 = ones(SIZE, 1);
    
    % iteration
    for iter = 1:N
        x1 = x2;
        for i = 1:SIZE
            s = 0;
            if i ~= 1 && i ~= SIZE
                low = i - 1;
                high = i + 1;
            else
                if i == 1
                    low = i;
                    high = i + 1;
                end
                if i == SIZE
                    low = i - 1;
                    high = i;
                end
            end
            for j = low:high
                if i ~= j
                    s = s + A(i, j) * x2(j);
                end
            end
            x2(i) = (omega * (b(i) - s) + (1 - omega) * A(i, i) * x2(i)) / A(i, i);
        end

        % compute error
        e = norm(x2 - exact, inf);
        error(iter, 1) = e;
        if e < tolerance
            break
        end
    end
end