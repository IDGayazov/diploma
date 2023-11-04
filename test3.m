clc; clear;

m = 10;
n = 15;
density = 0.2;
A = sprandsym(n, density, 0.8, 1);

% Проверка на положительную определенность
eigenvalues = eig(full(A)); % Получаем собственные значения матрицы
if all(eigenvalues > 0)
    disp('Матрица A положительно определена.');
else
    disp('Матрица A не положительно определена.');
end

B = A;

M = eye(n);

tau = 1;