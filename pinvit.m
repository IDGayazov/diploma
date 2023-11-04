function [X, lambda, iter] = pinvit(A, M, m, n, tau, eps, L, U, P)
%PINVIT Нахождение собственных значений и векторов матрицы
%  Предобусловленный метод для нахождения m собственных значений и векторов
% пары симметричных и положительно определенных матриц A, M
% lambda - вектор результатов собственных значений
% X - матрица, столбцы которой соответствуют собственным векторам
% iter - вектор состоящий из количества итерации для достижения
% точности eps

lambda = zeros(1, m);
X = zeros(n, m);
iter = zeros(1, m);

max_iter = 100000;
B = A; %

for s = 1:m
    % начальная итерация
    x = randn(n, 1);
    if s > 1
       x = ort(x, X, s - 1);
    end
    x = x / (sqrt(dot(M*x, x)));
    lam = dot(A*x, x);
    lam1 = lam - 2*eps;
    count = 0;
    while count ~= max_iter && abs(lam1 - lam) >= eps
        r = A*x - lam*M*x;
        u = B\r;
%         y = L \ (P*r);
%         u = U \ y;
        v = u;
        for j = 1:s-1
            Mu = M*u;
            alpha = dot(Mu, X(:, j));
            xp = X(:, j);
            v = v - alpha*xp;
        end
        x1 = x - tau*v;
        x2 = x1 / sqrt(dot(M*x1, x1));
        lam1 = lam;
        x = x2;
        lam = dot(A*x2, x2);
        count = count + 1;
    end
    iter(s) = count;
    lambda(s) = lam;
    X(:, s) = x;
end
end