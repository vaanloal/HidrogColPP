function [y] = dH_vap_col(T)

    % Válido de 120 < T < 600 oC
    % Dados obtidos através da equação de Vetere + Watson
    % T em degC
    % Vem em kJ/mol

    % Goodness of fit:
    %   SSE: 4.831
    %   R-square: 0.9999
    %   Adjusted R-square: 0.9999
    %   RMSE: 0.1582

    p1 = -5.643e-14;
    p2 = 1.104e-10;
    p3 = -8.731e-08;
    p4 = 3.531e-05;
    p5 = -0.007725;
    p6 = 0.7859;
    p7 = 66.36;

    y = p1 * T.^6 + p2 * T.^5 + p3 * T.^4 + p4 * T.^3 + p5 * T.^2 + p6 * T + p7;

end
