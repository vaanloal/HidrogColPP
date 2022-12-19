clc; clear all; close all; format long

% addpath("aux_fun path")
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');

%%%%%%%% MAIN INPUTS %%%%%%%%
treaction = 2.5; % h
toutros = 2; % h
Tmin = 140; % oC
Tmax = 200; % oC
TdegC0 = 178; % oC
Tfiltro_halt_degC = 150; % oC
Tcoolant_degC = mean([20, 50]); % oC
U_cooling = 400; % 500 - 1200 W/m2/K
tempo_transicao = 5; % min

tr_min = 0.5;
tr_max = 2.5;
%%%%%% END MAIN INPUTS %%%%%%

TK0 = TdegC0 + 273.15; % K
frac = [0.487772819 0.046512298 0.028606858 0.164908937 0.035763384 0.138086399 0.098349306 0];
frac = [0.309422684289359	0.0458880722748179	0.0625631193684445	0.286655191060494	0.0583190052412689	0.145108410333857	0.0895904633869377	0.00245305404482194];
% [frac] = getRandFrac();
% AA(1)   DIA(2) TEA(3) PAA(4)  DEA(5) NEA(6)  PIA(7) DPIAs(8)
% Proporções dos Ácidos Resínicos

spec.HeadSpace = 0.2; %
spec.TC = treaction + toutros; % h
spec.TLplus = 6.5; % h
spec.ProducaoAnual = 4785; % ton/ano
spec.SemanasAno = 48;
spec.HorasDia = 24; % h
spec.DiasSemana = 7;
spec.HD_racio = 1.3; % Razão Altura/Diâmetro do reator

[ReactorOut] = getVolumeReactor(frac, TdegC0, spec);

% Simulação Inicial
odeFun1 = @(t, x) rhscorr(t, x, TK0);
[t1, C1] = ode15s(odeFun1, [0, treaction * 60], ReactorOut.C0_sorted');

% Métricas da Simulação I
[metricas_simI] = MetricasReacao(C1);

% Quenching
Tcoolant = Tcoolant_degC + 273.15;
Tfiltro_halt = Tfiltro_halt_degC + 273.15;
A = ReactorOut.Area_Camisa; % m2
Cp_mix = 2.145114739 * 1000; % J/(kg oC)
m_mix = ReactorOut.MassaMix; % kg

par(1) = Tcoolant;
par(2) = m_mix;
par(3) = Cp_mix;
par(4) = U_cooling;
par(5) = A;
par(6) = tempo_transicao;
par(7) = Tfiltro_halt;

C0_quench = C1(end, :);
T0_quench = TK0;
x0_quench = [C0_quench, T0_quench]';

options = odeset('Events', @(t, x) quenchingEvent(t, x, Tfiltro_halt));
[t2, C2] = ode15s(@(t, x) rhs_quenching(t, x, par), [0, 5000], x0_quench, options);

t_all = [t1; t2(2:end) + t1(end)];
C_all = [C1; C2(2:end, 1:end - 1)];
T_all = [ones(size(t1, 1), 1) * TK0; C2(2:end, end)] - 273.15;

% Métricas da Simulação II
[metricas_simII] = MetricasReacao(C_all);

%% PLOT KINETICS
nomes = ["$C_{AA}$", "$C_{DIA}$", "$C_{TEA}$", "$C_{PAA}$", "$C_{DEA}$", "$C_{NEA}$", "$C_{PIAs}$", "$C_{DPIAs}$"];
all_marks = {'o', '+', '*', 'h', 'x', 's', 'd', '^', 'v', '>', '<', 'p', '.'};
hfig = figure("Color", "w", Position = [4.868333333333334, 10.583333333333336, 20.002500000000005, 13.017500000000002]);
tiledlayout(1, 2)
lin_width = 1;
nexttile
hold on
plot(t_all / 60, C_all(:, 1), LineWidth = lin_width, DisplayName = "$C_{AA}$", Marker = all_marks{1}, MarkerIndices = [1 30 40 50 60 70 80 90])
plot(t_all / 60, C_all(:, 2), LineWidth = lin_width, DisplayName = "$C_{DIA}$", Marker = all_marks{2}, MarkerIndices = [1 30 40 50 60 70 80 90])
plot(t_all / 60, C_all(:, 3), LineWidth = lin_width, DisplayName = "$C_{TEA}$", Marker = all_marks{3}, MarkerIndices = [1 30 40 50 60 70 80 90])
plot(t_all / 60, C_all(:, 4), LineWidth = lin_width, DisplayName = "$C_{PAA}$", Marker = all_marks{4}, MarkerIndices = [1 30 40 50 60 70 80 90])
plot(t_all / 60, C_all(:, 5), LineWidth = lin_width, DisplayName = "$C_{DEA}$", Marker = all_marks{5}, MarkerIndices = [1 30 40 50 60 70 80 90])
plot(t_all / 60, C_all(:, 6), LineWidth = lin_width, DisplayName = "$C_{NEA}$", Marker = all_marks{6}, MarkerIndices = [1 30 40 50 60 70 80 90])
plot(t_all / 60, C_all(:, 7), LineWidth = lin_width, DisplayName = "$C_{PIAs}$", Marker = all_marks{7}, MarkerIndices = [1 30 40 50 60 70 80 90])
plot(t_all / 60, C_all(:, 8), LineWidth = lin_width, DisplayName = "$C_{DPIAs}$", Marker = all_marks{8}, MarkerIndices = [1 30 40 50 60 70 80 90])
legend(NumColumns = 2, Interpreter = "latex", Location = "best", Box = "off")
xline(treaction + tempo_transicao / 60, LineStyle = "--", Color = [0.15 0.15 0.15])
text((treaction + tempo_transicao / 60) * 0.975, 1 * 0.5, 'Despressurização', Rotation = 90, FontName = 'Georgia')
set(gca, "Box", "on")
xlabel("$t\; \rm[h]$", Interpreter = "latex")
ylabel("$C_i\;\rm[mol/L]$", Interpreter = "latex")

ylim([0, inf])
xlim([-0.01, 1.05 * max(t_all / 60)])

nexttile
plot(t_all / 60, T_all, LineWidth = 2); hold on
xline(treaction + tempo_transicao / 60, LineStyle = "--", Color = [0.15 0.15 0.15])
text((treaction + tempo_transicao / 60) * 0.975, 1 * min(T_all), 'Despressurização', Rotation = 90, FontName = 'Georgia')
xlabel("$t\; \rm[h]$", Interpreter = "latex")
ylabel("$T\; \rm[^oC]$", Interpreter = "latex")
set(gca, "FontSize", 12)

SetFigureDefaults2(hfig, 0.65, 20, 12, gca)
%%

figure(Color = "w")
tiledlayout("flow")
nexttile
plot(t_all / 60, C_all, LineWidth = 2); hold on
xline(treaction + tempo_transicao / 60, LineStyle = "--", Color = [0.15 0.15 0.15])
text((treaction + tempo_transicao / 60) * 0.975, 0.6 * max(C_all, [], 'all'), 'Despressurização', Rotation = 90, FontName = 'Georgia')
ylim([0, inf])
legend(nomes, Interpreter = "latex", NumColumns = 2, Location = "best", Box = "off")
xlabel("$t\; [\rm h]$", Interpreter = "latex")
ylabel("$C_i\;\rm[mol/L]$", Interpreter = "latex")
set(gca, "FontSize", 12)

nexttile
plot(t_all / 60, T_all, LineWidth = 2); hold on
xline(treaction + tempo_transicao / 60, LineStyle = "--", Color = [0.15 0.15 0.15])
text((treaction + tempo_transicao / 60) * 0.975, 1 * min(T_all), 'Despressurização', Rotation = 90, FontName = 'Georgia')
xlabel("$t\; \rm[h]$", Interpreter = "latex")
ylabel("$T\; \rm[^oC]$", Interpreter = "latex")
set(gca, "FontSize", 12)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Otimizar a temperatura de reação para um dado tempo estipulado %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

format short
tt = linspace(tr_min, tr_max, 30);
Tstar = ones(length(tt), 1);
obj_iso = ones(length(tt), 1);
Pureza_iso = ones(length(tt), 1);

tic

for i = 1:length(tt)
    treaction = tt(i); % h
    fun = @(T) OptimizeReactorIso(treaction, frac, spec, T);
    [Tstar(i), obj_iso(i)] = fminbnd(fun, Tmin, Tmax);
    metricas_Tstar = OptimizeReactorIso2(treaction, frac, spec, Tstar(i));
    Pureza_iso(i) = metricas_Tstar.Pureza;
end

time_optimization = toc;

figure(Color = "w")
yyaxis left
plot(tt, -1 * obj_iso * 100, LineWidth = 2, DisplayName = "OBJ"); hold on
plot(tt, Pureza_iso * 100, LineWidth = 2, DisplayName = "Pureza (%)")
yline(92, DisplayName = "Referência 1", LineStyle = "--")
yline(90, DisplayName = "Referência 1", LineStyle = ":")

legend(Box = "off")
ylabel("OBJ (%)")
xlabel("t [h]")

yyaxis right
plot(tt, Tstar, LineWidth = 2, DisplayName = "Temperatura")
ylabel("Temperatura")
xlabel("t [h]")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Adicionar Incerteza ao passo anterior com Monte Carlo %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

no_runs = 25;
t_MC = linspace(tr_min, 8, 30);

TstarMC = ones(length(t_MC), no_runs);
obj_isoMC = ones(length(t_MC), no_runs);
Pureza_isoMC = ones(length(t_MC), no_runs);
frac_store = zeros(no_runs, length(frac));

tic

for j = 1:no_runs
    % Monte Carlo Loop

    [frac] = getRandFrac();
    frac_store(j, :) = frac;

    for i = 1:length(t_MC)
        treaction = t_MC(i); % h
        fun = @(T) OptimizeReactorIso(treaction, frac, spec, T);
        [TstarMC(i, j), obj_isoMC(i, j)] = fminbnd(fun, Tmin, Tmax);
        metricas_TstarMC = OptimizeReactorIso2(treaction, frac, spec, TstarMC(i, j));
        Pureza_isoMC(i, j) = metricas_TstarMC.Pureza;
    end

    j
end

time_MonteCarlo = toc;

T_MC_mean = mean(TstarMC, 2);
T_MC_std = std(TstarMC, [], 2);

OBJ_MC_mean = -100 * mean(obj_isoMC, 2);
OBJ_MC_std = -100 * std(obj_isoMC, [], 2);

Pureza_MC_mean = 100 * mean(Pureza_isoMC, 2);
Pureza_MC_std = 100 * std(Pureza_isoMC, [], 2);

%% Um único tempo (Distribuição de Temperaturas isotérmicas)

no_runs = 1000;
t_MC = 2.5; % h

TstarMC = ones(length(t_MC), no_runs);
obj_isoMC = ones(length(t_MC), no_runs);
Pureza_isoMC = ones(length(t_MC), no_runs);
frac_store = zeros(no_runs, length(frac));

tic

for j = 1:no_runs
    % Monte Carlo Loop

    [frac] = getRandFrac();
    frac_store(j, :) = frac;

    treaction = t_MC; % h
    fun = @(T) OptimizeReactorIso(treaction, frac, spec, T);
    [TstarMC(j), obj_isoMC(j)] = fminbnd(fun, Tmin, Tmax);
    metricas_TstarMC = OptimizeReactorIso2(treaction, frac, spec, TstarMC(j));
    Pureza_isoMC(j) = metricas_TstarMC.Pureza;
    j
end

time_MonteCarlo = toc;

T_MC_mean = mean(TstarMC, 2);
T_MC_std = std(TstarMC, [], 2);

OBJ_MC_mean = -100 * mean(obj_isoMC, 2);
OBJ_MC_std = -100 * std(obj_isoMC, [], 2);

Pureza_MC_mean = 100 * mean(Pureza_isoMC, 2);
Pureza_MC_std = 100 * std(Pureza_isoMC, [], 2);

%%

face_alpha = 1;
edge_alpha = 0.3;
edge_color = [0 0 0]; %[0.96 0.96 0.96];
face_color = [0.90, 0.90, 0.90];
linwidth = 0.3;

methods = {'fd'; 'scott'; 'sturges'};
    idx = 3;
    n = calcnbins(TstarMC, methods{idx});
    figure(Color = "W", Position = [422, 488, 1109, 420])
    tiledlayout(1, 2, 'TileSpacing', 'Compact')
    nexttile
    h1 = histfit(TstarMC, n, 'kernel', FaceColor = face_color);
    set(h1(1), 'FaceAlpha', face_alpha);
    set(h1(1), 'FaceColor', face_color)
    set(h1(1), 'EdgeColor', edge_color);
    set(h1(1), 'LineWidth', linwidth);
    xlabel("$T\; \rm[^oC]$", Interpreter = "latex")
    ylabel("$\rm Frequ\hat{e}ncia$", Interpreter = "latex")
    set(gca, "FontSize", 12)
    nexttile

    n = calcnbins(Pureza_isoMC, methods{idx});
    h1 = histfit(Pureza_isoMC * 100, n, 'kernel', FaceColor = face_color);
    set(h1(1), 'FaceAlpha', face_alpha);
    set(h1(1), 'FaceColor', face_color)
    set(h1(1), 'EdgeColor', edge_color);
    set(h1(1), 'LineWidth', linwidth);
    xlabel('$\chi\;\rm [\%]$', Interpreter = 'latex')
    ylabel("$\rm Frequ\hat{e}ncia$", Interpreter = "latex")
    set(gca, "FontSize", 12)

    exportgraphics(gcf, 'pureza_tempo_dist_iso.pdf')

    %%

    %%
    figure(Color = "w")
    yyaxis left
    plot(t_MC, OBJ_MC_mean, LineWidth = 2, DisplayName = "\psi", LineStyle = "-"); hold on
    plot(t_MC, OBJ_MC_mean - OBJ_MC_std, LineWidth = 1, DisplayName = "\psi - \sigma", LineStyle = ":"); hold on
    plot(t_MC, OBJ_MC_mean + OBJ_MC_std, LineWidth = 1, DisplayName = "\psi + \sigma", LineStyle = ":"); hold on

    plot(t_MC, Pureza_MC_mean, LineWidth = 2, DisplayName = "\eta (%)", LineStyle = "--"); hold on
    plot(t_MC, Pureza_MC_mean - Pureza_MC_std, LineWidth = 1, DisplayName = "\eta (%) - \sigma", LineStyle = ":", Marker = "none"); hold on
    plot(t_MC, Pureza_MC_mean + Pureza_MC_std, LineWidth = 1, DisplayName = "\eta (%) + \sigma", LineStyle = ":", Marker = 'none'); hold on

    yline(92, DisplayName = "Referência 1", LineStyle = "--")
    yline(90, DisplayName = "Referência 2", LineStyle = ":")

    legend(Box = "off")
    ylabel("\psi, \eta (%)")
    xlabel("t [h]")

    yyaxis right
    plot(t_MC, T_MC_mean, LineWidth = 2, DisplayName = "T", LineStyle = "-."); hold on
    plot(t_MC, T_MC_mean - T_MC_std, LineWidth = 1, DisplayName = "T - \sigma", LineStyle = ":"); hold on
    plot(t_MC, T_MC_mean + T_MC_std, LineWidth = 1, DisplayName = "T + \sigma", LineStyle = ":"); hold on
    ylabel("Temperatura [^oC]")
    xlabel("t [h]")

    set(gca, "FontSize", 12)