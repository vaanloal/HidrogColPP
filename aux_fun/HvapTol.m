function [Y] = HvapTol(T_K)
    Tc = 318.60 + 273.15; % K

    T_r = T_K / Tc;
    A = 54643000;
    B = 0.76764;
    C = -0.62056;
    D = 0.25935;
    E = 0;
    F = 0;
    G = 0;

    Y = A * (1 - T_r).^(B + C * T_r + D * T_r.^2 + E * T_r.^3); % J/kmol
    Y = Y / 1000; %kJ/kmol
end
