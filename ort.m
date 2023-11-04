function res = ort(x, X, s)
% построение ортогонального вектора (ортогонолизация Грама-Шмидта)
res = x;
if s ~= 0
    for i = 1:s
        bb = dot(X(:, i), X(:, i));
        if  bb ~= 0
            res = res - dot(X(:, i), x) * X(:, i) / bb;
        end
    end
end
end