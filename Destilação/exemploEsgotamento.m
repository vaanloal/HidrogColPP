clc; clear all; close all
%% INPUTS

P = 15; %mmHg
per_vac = (P / 760) * 100;

F = 8228.919592; % mol/h
xF = 0.739860073; % xileno
xD = 0.999;
xB = 0.0028;

q = 1.0001;
eff = 1;
murphree_eff = 0.85;
Lp = 0.5; % [m] 0.3048 * 2; % Espaçamento entre pratos entre 0.15 e 1 metros

Mw_col = 304.5;
Mw_xil = 106.2;
dH_vap_col_fix = 73.914; % kJ/mol
Cp_col = 2.03059800; %  kJ/(kg*K)

addpath("D:\Engenharia Química\5º Ano\Projeto de Processo\Codigos\propriedades\aux_fun")

Ant.Col = [7.06195 1850.91 66.1306]; % [mmHg] válida de 173.5 oC a 558.85
Ant.Xil = [7.15471 1553.95 225.23]; % [mmHg] válida de 13.23 oC a 343.11 oC
P_vap_Col = @(TdegC) 10 .^ (Ant.Col(1) - Ant.Col(2) ./ (TdegC + Ant.Col(3)));

P_vap_Xyl = @(TdegC) 10 .^ (Ant.Xil(1) - Ant.Xil(2) ./ (TdegC + Ant.Xil(3)));
T_vap_Col = @(P) (Ant.Col(2) - (Ant.Col(1) - log10(P)) * Ant.Col(3)) / (Ant.Col(1) - log10(P));
T_vap_Xyl = @(P) (Ant.Xil(2) - (Ant.Xil(1) - log10(P)) * Ant.Xil(3)) / (Ant.Xil(1) - log10(P));

[alpha_mean, alpha] = GetAlpha(Ant, P, 0);
% exportgraphics(gca, "volatilidade_plot.pdf")

alpha_bar = min(alpha);

%% Coluna de Esgotamento
% Curva de equilibrio

figure(Color = "w", Position = [573, 501.7, 590.7, 462.3])
plot([0, 1], [0, 1], "k-", LineWidth = 1)
xlabel("$x$", Interpreter = "latex")
ylabel("$y$", Interpreter = "latex")
set(gca, "XTick", 0:0.1:1)
set(gca, "YTick", 0:0.1:1)
grid on
ylim([0, 1])
xlim([0, 1])
hold on

% Curva de equilibrio
y_eq = @(x) alpha_bar * x ./ (1 - (1 - alpha_bar) .* x);
fplot(y_eq, [0, 1], Color = "r", LineWidth = 1)

% Linha dos q's
q_line = @(x) q * x ./ (q - 1) - xF / (q - 1);

% ponto de interseção com a linha dos q's
if q == 1 % evitar valor infinito (divisão por 0)
    x_pinch_point = xF;
else
    x_pinch_point = fsolve(@(x) q_line(x) - y_eq(x), 0.5, optimset('Display', 'off'));
end

y_pinch_point = y_eq(x_pinch_point);

if q == 1 % evitar valor infinito (divisão por 0)
    x_intersect = xF;
else
    x_intersect = fsolve(@(x) q_line(x) - xD, xF, optimset('Display', 'off'));
end

y_intersect = xD;

% q-line
if q > 1
    fplot(q_line, [xF, x_intersect], Color = "m", LineWidth = 1)

elseif q == 1
    plot([xF, xF], [xF, x_intersect], Color = "m", LineWidth = 1)

elseif q > 0 && q < 1
    fplot(q_line, [0, xF], Color = "m", LineWidth = 1)

elseif q == 0
    plot([xF, 0], [xF, xF], Color = "m", LineWidth = 1)

else
    fplot(q_line, [0, xF], Color = "m", LineWidth = 1)
end

% Top Line
plot([x_intersect, xD], [xD, xD], Color = "g", LineWidth = 1)

% Bot Line
plot([xB, x_intersect], [xB, y_intersect], Color = "g", LineWidth = 1)

% Desenhar os estágios
x_eq = @(y) y ./ (alpha_bar + (1 - alpha_bar) .* y);

BottomSlope = (y_intersect - xB) / (x_intersect - xB);
BottomInter = (BottomSlope - 1) * xB;
x_bottom = @(y) (y - BottomInter) / BottomSlope;

p = polyfit([xB, x_intersect], [xB, y_intersect], 1);
x_bottom2 = @(y) (y + p(2)) ./ p(1);
y_bottom = @(x) p(1) * x + p(2);

x_top_1 = x_intersect;
y_top_1 = y_intersect;

loop = 0;

while x_top_1 > xB && loop < 200

    y_top_2 = y_top_1;
    x_top_2 = x_eq(y_top_2);
    x_top_2 = murphree_eff * (x_top_2 - x_top_1) + x_top_1;
    x_top_3 = x_top_2;
    y_top_3 = y_bottom(x_top_3);

    plot([x_top_1 x_top_2], [y_top_1 y_top_2], "b", LineWidth = 1)
    plot([x_top_2 x_top_3], [y_top_2 y_top_3], "b", LineWidth = 1)

    x_plot_top = x_top_1;
    x_top_1 = x_top_3;
    y_top_1 = y_top_3;

    loop = loop + 1;

    if loop >= 199
        warning("Número de estágios não convergiu!")
    end

end

set(gca, 'FontSize', 12)

if y_eq(x_intersect) < xD
    warning("Separação inviável. Especificar novas frações no destilado.")
    beep
    close all
else
    NoStages = loop - (x_top_2 - xB) / (x_top_2 - x_plot_top);
end

%% Calcular Calores e Temperaturas de Topo e da Base +
B =- ((F * xD - F * xF) / (xB - xD)); % caudal do resíduo [mol/h]
D = ((F * xB - F * xF) / (xB - xD)); % caudal do destilado [mol/h]
L = D * 0; % Refluxo [mol/h]
V = L + D; % Caudal de Vapor above feed  [mol/h]
Vbar = V - F + F * q; % Caudal de Vapor bellow feed [mol/h]

funDew = @(T) (P * xD ./ P_vap_Xyl(T) + P * (1 - xD) ./ P_vap_Col(T)) - 1; % dew point
funBub = @(T) (xB .* P_vap_Xyl(T) / P + (1 - xB) .* P_vap_Col(T) / P) - 1; % bubble point
T_topo = fsolve(funDew, 150, options);
T_base = fsolve(funBub, 180, options);

dH_vap_xil_topo = dH_vap_xil(T_topo + 273.15); % kJ/mol
dH_vap_xil_base = dH_vap_xil(T_base + 273.15); % kJ/mol

dH_vap_col_topo = dH_vap_col(T_topo); % kJ/mol
dH_vap_col_base = dH_vap_col(T_base); % kJ/mol

% Para o caso de condensador total
Qc = V * (dH_vap_xil_topo * xD + (1 - xD) * dH_vap_col_topo) / 3600; % [kJ/s ==> kW] REMOVER
Qr = Vbar * (dH_vap_xil_base * xB + (1 - xB) * dH_vap_col_base) / 3600; % [kJ/s ==> kW] ADICIONAR

%% Altura da coluna
% cada prato tem cerca de 2 ft (2*0.3048 m) de altura (Douglas)
margem = 0.15; % adicionar esta percentagem para o topo e base da coluna
H = ((ceil(NoStages) - 1) * Lp) * (1 + margem); % [m]

%% Diâmetro da coluna
RMM_feed = xF * Mw_xil + (1 - xF) * Mw_col; % [g/mol]
RMM_top = xD * Mw_xil + (1 - xD) * Mw_col; % [g/mol]
RMM_bot = xB * Mw_xil + (1 - xB) * Mw_col; % [g/mol]

mass_flow_feed = (xF * F * Mw_xil + (1 - xF) * F * Mw_col) / 1000; % [kg/h]
mass_flow_top = (xD * D * Mw_xil + (1 - xD) * D * Mw_col) / 1000; % [kg/h]
mass_flow_bot = (xB * B * Mw_xil + (1 - xB) * B * Mw_col) / 1000; % [kg/h]

%% Topo da coluna
rho_L_topo = xD * rho_xyl(T_topo) + (1 - xD) * rho_col(T_topo); % [kg/m^3] média entre a base e o topo
rho_v_topo = P * 133.322368 * RMM_top / (8.3145 * T_base * 1000); % [kg/m^3] média entre a base e o topo

% rho do gás e do liquido
rho_L_base = xB * rho_xyl(T_base) + (1 - xB) * rho_col(T_base); % [kg/m^3] média entre a base e o topo
rho_v_base = P * 133.322368 * RMM_bot / (8.3145 * T_base * 1000); % [kg/m^3] média entre a base e o topo

rho_L = mean([rho_L_base rho_L_topo]);
rho_v = mean([rho_v_base rho_v_topo]);

% correlação de Souders Brown
uv = (-0.171 * Lp ^ 2 + 0.27 * Lp - 0.047) * ((rho_L - rho_v) / rho_v) ^ 0.5; % velocidade de vapor máxima [m/s]

V_max = max([V * RMM_top Vbar * RMM_bot]); % [g/h]
V_w = V_max / 3600/1000; % [kg/s]

% Diâmetro da coluna
Dc = (4 * V_w / (pi * rho_v * uv)) ^ 0.5; % [m]
A_active = (pi / 4) * Dc ^ 2;

% Adicionar 12% de downcomers
Ac = A_active * 1.12; % [m2]
Dc = (4 * Ac / pi) ^ 0.5;

%% Custos da Coluna
horas_operacao = 8064;
anos = 3;
CEPCI2022 = 832.6;
CEPCI2017 = 567.5;
MES2022 = 2171.6 * 1.1 * 1.05;
CustoEletricidade = 0.1572; %€/kWh
CustoEletricidadeFixo = 0.7780; %€/dia
PriceCoolingWater = (CEPCI2022 / CEPCI2017) * 0.378/1000000; % $/kJ

Fm = 3.67; % stainless steal
[CustoColuna, InstalacaoColuna] = TowerCost(Dc, H, P, Fm, MES2022);
Ft = 0;
Fm = 1.7; % stainless steal
CustoPratos = TraysCost(Dc, H, Fm, Ft, Lp, MES2022);

%% Outputs
exportgraphics(gca, "mccabe_thiele.pdf", resolution = 600)

fprintf('\n______________________________________________________')
fprintf('\n                        OUTPUTS                     \n')
fprintf('______________________________________________________\n')
fprintf('(1)  Número de andares  . . . . . . . . .  %.2f [-]\n', ceil(NoStages))
fprintf('(2)  Altura da coluna  . . . . . . . . . . %.2f [m]\n', H)
fprintf('(3)  Espaçamento dos pratos (Lp)  . . . .  %.2f [m]\n', Lp)
fprintf('(4)  Diâmetro da coluna  . . . . . . . . . %.2f [m]\n', Dc)
fprintf('(5)  Caudal de destilado (D)  . . . . . .  %.2f [mol/h]\n', D)
fprintf('(6)  Caudal do resíduo (B)  . . . . . . .  %.2f [mol/h]\n', B)
fprintf('(7)  Caudal de destilado (D_kg)  . . . . . %.2f [kg/h]\n', mass_flow_top)
fprintf('(8)  Caudal do resíduo (B_kg)  . . . . . . %.2f [kg/h]\n', mass_flow_bot)
fprintf("(9)  Temperatura do destilado (TD)  . . .  %.2f [oC]\n", T_topo)
fprintf("(10) Temperatura do resíduo (TB)  . . . .  %.2f [oC]\n", T_base)
fprintf("(11) Calor do revaporizador (Qr)  . . . .  %.2f [kW]\n", Qr)
fprintf("(13) Calor do condensador (Qc)  . . . . .  %.2f [kW]\n", Qc)
fprintf('______________________________________________________\n')

REC_TOP_XYL = (D * xD) / (F * xF);
REC_TOP_COL = (D * (1 - xD)) / (F * (1 - xF));
vec2excel = [REC_TOP_XYL REC_TOP_COL]
