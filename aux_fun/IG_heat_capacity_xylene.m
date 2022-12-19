function [Y] = IG_heat_capacity_xylene(TdegC)

% Entalpia do gás ideal em J/kmol·K


if TdegC > 1500-273.15
    warning('A propriedade está a ser extrapolada')
elseif TdegC < 298.15-273.15
    warning('A propriedade está a ser extrapolada')
end

T = TdegC + 273.15;

A = 92948;
B = 242050;
C = 800.71;
D = 113980;
E = 2312;
F = 0;
G = 0;

Y = A + B * ((C ./ T) ./ (sinh(C ./ T))).^2 +...
    D * ((E ./ T) ./ (cosh(E ./ T))).^2;

end


