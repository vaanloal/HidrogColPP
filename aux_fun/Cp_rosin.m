function [cp] = Cp_rosin(T)
% Dados experimentais
% Válida de 60 oC a 260 oC
% T = oC

Tmax = 260;
Tmin = 60;

idxmax = find(T >= Tmax);
T(idxmax) = Tmax;

idxmin = find(T <= Tmin);
T(idxmin) = Tmin;

if ~isempty(idxmin)
    warning('A propriedade está a ser extrapolada para baixo!!')
elseif ~isempty(idxmax)
      warning('A propriedade está a ser extrapolada para cima!!')
end



cp = 1.775226206141090 + 0.001696343374692 * T; % kJ/(kg.oC)

end

