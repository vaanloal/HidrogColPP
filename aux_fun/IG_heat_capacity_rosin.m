function y = IG_heat_capacity_rosin(TdegC)
% J/(mol.K)
% Válida de 150 a 600 oC
% Valores obtidos através do método de Joback

if TdegC > 600
    warning('A propriedade está a ser extrapolada')
elseif TdegC < 150
    warning('A propriedade está a ser extrapolada')
end

y = - 6.363e-07.*TdegC.^2 + 0.001069.*TdegC + 0.433;
y = y*1000;

end