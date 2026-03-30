function y = u0(x1, x2)
    condition = ((x1.^2 + x2.^2) <= 0.64 & (x1.^2 + x2.^2) >= 0.04);
    y = zeros(size(x1));
    y(condition) = 10;
end
