function [Y] = HvapBenz(T_K)

    Tc = 288.90 + 273.15; % K
    T_r = T_K / Tc;

    A = 50007000;
    B = 0.65393;
    C = -0.27698;
    D = 0.029569;
    E = 0;
    F = 0;
    G = 0;

    Y = A * (1 - T_r).^(B + C * T_r + D * T_r.^2 + E * T_r.^3); % J/kmol
    Y = Y / 1000; %kJ/kmol
end
