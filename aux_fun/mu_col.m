function [mu_col] = mu_col(T)
A = -28.9145;
B = 1987.6175;
C = 4.0051;

mu_col = exp(A + B./ T + C * log(T)); % cP
end