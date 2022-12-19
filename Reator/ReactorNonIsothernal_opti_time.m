clc; clear all; close all; format long
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');

% Adicionar o CasADi ao path do Matlab
% addpath("CasADi path")
% addpath("aux_fun path")

%%%%%%%% MAIN INPUTS %%%%%%%%
treaction = 2.5; % h
toutros = 2; % h
Tmin = 140; % oC
Tmax = 200; % oC
TdegC0 = 175; % oC
Tfiltro_halt_degC = 150; % oC
Tcoolant_degC = mean([20, 50]); % oC
U_cooling = 400; % 500 - 1200 W/m2/K
tempo_transicao = 5; % min

tr_min = 0.5;
tr_max = 3.5;

spec.HeadSpace = 0.2; %
spec.TC = treaction + toutros; % h
spec.TLplus = 6.5; % h
spec.ProducaoAnual = 4785; % ton/ano
spec.SemanasAno = 48;
spec.HorasDia = 24; % h
spec.DiasSemana = 7;
spec.HD_racio = 1.3; % Razão Altura/Diâmetro do reator

N = 100; % Número de Controlos
RendimentoGOAL = 0.94;

%%%%%% END MAIN INPUTS %%%%%%
%_________________________________________________________________________

% Proporções dos Ácidos Resínicos
% AA(1)   DIA(2) TEA(3) PAA(4)  DEA(5) NEA(6)  PIA(7) DPIAs(8)
frac = [0.487772819 0.046512298 0.028606858 ...
        0.164908937 0.035763384 0.138086399 0.098349306 0];

% Obter Condições Iniciais
[ReactorOut] = getVolumeReactor(frac, TdegC0, spec);
[IC] = getIC(treaction, TdegC0, ReactorOut, N);
[metricas_sim] = MetricasReacao(IC);

%% Otimização Dinâmica + Monte Carlo min(Time)
n_runs = 500;

for j = 1:n_runs
    [frac] = getRandFrac();

    [metricas_opti(j), res(j)] = RendDynTime(frac, TdegC0, spec, N, Tmin, ...
        Tmax, tr_min, tr_max, IC, RendimentoGOAL);

end

%% Plots

all_marks = {'o', '+', '*', 'h', 'x', 's', 'd', '^', 'v', '>', '<', 'p', '.'};
hfig = figure("Color", "w", Position = [4.868333333333334, 10.583333333333336, 20.002500000000005, 13.017500000000002]);
lin_width = 1;
hold on
plot(res.tgrid, res.states(:, 1), LineWidth = lin_width, DisplayName = "$C_{AA}$", Marker = all_marks{1}, MarkerIndices = (N / 10:N / 10:N))
plot(res.tgrid, res.states(:, 2), LineWidth = lin_width, DisplayName = "$C_{DIA}$", Marker = all_marks{2}, MarkerIndices = (N / 10:N / 10:N))
plot(res.tgrid, res.states(:, 3), LineWidth = lin_width, DisplayName = "$C_{TEA}$", Marker = all_marks{3}, MarkerIndices = (N / 10:N / 10:N))
plot(res.tgrid, res.states(:, 4), LineWidth = lin_width, DisplayName = "$C_{PAA}$", Marker = all_marks{4}, MarkerIndices = (N / 10:N / 10:N))
plot(res.tgrid, res.states(:, 5), LineWidth = lin_width, DisplayName = "$C_{DEA}$", Marker = all_marks{5}, MarkerIndices = (N / 10:N / 10:N))
plot(res.tgrid, res.states(:, 6), LineWidth = lin_width, DisplayName = "$C_{NEA}$", Marker = all_marks{6}, MarkerIndices = (N / 10:N / 10:N))
plot(res.tgrid, res.states(:, 7), LineWidth = lin_width, DisplayName = "$C_{PIAs}$", Marker = all_marks{7}, MarkerIndices = (N / 10:N / 10:N))
plot(res.tgrid, res.states(:, 8), LineWidth = lin_width, DisplayName = "$C_{DPIAs}$", Marker = all_marks{8}, MarkerIndices = (N / 10:N / 10:N))
legend(NumColumns = 4, Interpreter = "latex", Location = "northoutside", Box = "off")
set(gca, "Box", "on")
xlabel("$t\; [h]$", Interpreter = "latex")
ylabel("$C_i\;\rm[mol/L]$", Interpreter = "latex")
SetFigureDefaults2(hfig, 0.75, 20, 12, gca)
ylim([0, inf])
xlim([-0.01, max(res.tgrid)])

figure(Color = "w")
stairs(res.tgrid(1:end - 1), res.controls(:, 1), "r-", LineWidth = lin_width, DisplayName = "$T$")
xlabel("$t\; [h]$", Interpreter = "latex")
ylabel("$T\; \rm[^oC]$", Interpreter = "latex")
SetFigureDefaults(gca)
axis padded

for i = 1:n_runs
    T_result(:, i) = res(i).controls;
    tgrid(:, i) = res(i).tgrid(1:end - 1);
    % plot(tgrid, T_result)
end

%%
figure(Color = "W")

for i = 1:n_runs
    plot(tgrid(:, i), T_result(:, i)); hold on
end

% plot(mean(tgrid,2), mean(T_result,2), Color=[0,0,0], LineWidth=2)
xlabel("$t\; [h]$", Interpreter = "latex")
ylabel("$T\; \rm[^oC]$", Interpreter = "latex")
axis padded
set(gca, "FontSize", 12)
exportgraphics(gca, 'Temp_profiles_min_time.pdf')

for i = 1:n_runs
    time_dist(i, 1) = res(i).tf;
    pur_dist(i, 1) = metricas_opti(i).Pureza;
end

%%
close all
face_alpha = 1;
edge_alpha = 0.3;
edge_color = [0 0 0]; %[0.96 0.96 0.96];
face_color = [0.90, 0.90, 0.90];
linwidth = 0.3;

methods = {'fd'; 'scott'; 'sturges'};
    idx = 3;
    n = calcnbins(time_dist, methods{idx});
    figure(Color = "W", Position = [422, 488, 1109, 420])
    tiledlayout(1, 2, 'TileSpacing', 'Compact')
    nexttile
    h1 = histfit(time_dist, n, 'kernel', FaceColor = face_color);
    set(h1(1), 'FaceAlpha', face_alpha);
    set(h1(1), 'FaceColor', face_color)
    set(h1(1), 'EdgeColor', edge_color);
    set(h1(1), 'LineWidth', linwidth);
    xlabel('$t\;\rm [h]$', Interpreter = 'latex')
    ylabel("$\rm Frequ\hat{e}ncia$", Interpreter = "latex")
    set(gca, "FontSize", 12)
    nexttile

    n = calcnbins(pur_dist, methods{idx});
    h1 = histfit(pur_dist * 100, n, 'kernel', FaceColor = face_color);
    set(h1(1), 'FaceAlpha', face_alpha);
    set(h1(1), 'FaceColor', face_color)
    set(h1(1), 'EdgeColor', edge_color);
    set(h1(1), 'LineWidth', linwidth);
    xlabel('$\chi\;\rm [\%]$', Interpreter = 'latex')
    ylabel("$\rm Frequ\hat{e}ncia$", Interpreter = "latex")
    set(gca, "FontSize", 12)

    exportgraphics(gcf, 'pureza_tempo_dist.pdf')
