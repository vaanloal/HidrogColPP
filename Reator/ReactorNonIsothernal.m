clc; clear all; close all; format long

% Adicionar o CasADi ao path do Matlab
% addpath("CasADi path")
% addpath("aux_fun path")

set(groot, 'defaultAxesTickLabelInterpreter', 'latex');

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
%%%%%% END MAIN INPUTS %%%%%%

TK0 = TdegC0 + 273.15; % K
frac = [0.487772819 0.046512298 0.028606858 0.164908937 0.035763384 0.138086399 0.098349306 0];
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

% tfinal = (m_mix * Cp_mix / (U_cooling * ReactorOut.Area_Camisa)) * log((150 - Tcoolant_degC)./(Tfiltro_halt_degC - Tcoolant_degC))

% ESTIMAR CAUDAL DE ÁGUA NECESSÁRIO
% Q_cooling = m_mix * Cp_mix * (T0_quench - Tfiltro)/1000; % kJ
% T_agua_out = 50; % oC
% T_agua_in  = 20; % oC
% Cp_water = 4.186; % kJ/(kg.oC)
% m_water = Q_cooling / (Cp_water * (T_agua_out - T_agua_in)); % kg

t_all = [t1; t2(2:end) + t1(end)];
C_all = [C1; C2(2:end, 1:end - 1)];
T_all = [ones(size(t1, 1), 1) * TK0; C2(2:end, end)] - 273.15;

% Métricas da Simulação II
[metricas_simII] = MetricasReacao(C_all);

nomes = ["$C_{AA}$", "$C_{DIA}$", "$C_{TEA}$", "$C_{PAA}$", "$C_{DEA}$", "$C_{NEA}$", "$C_{PIAs}$", "$C_{DPIAs}$"];

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
text((treaction + tempo_transicao / 60) * 0.975, 1.2 * min(T_all), 'Despressurização', Rotation = 90, FontName = 'Georgia')
xlabel("$t\; \rm[h]$", Interpreter = "latex")
ylabel("$T\; \rm[^oC]$", Interpreter = "latex")
set(gca, "FontSize", 12)

%%%%%%%%%%%%%%%%%%%%%%%%%
%% Otimização Dinâmica %%
%%%%%%%%%%%%%%%%%%%%%%%%%
% Econtrar o rendimento 'ótimo' para um dado tempo com Monte Carlo
%___________________________________________________________________%
%%
N = 100;

no_runs = 350;
t_MC = linspace(tr_min, 4, 30);

TstarMC = ones(length(t_MC), no_runs);
obj_isoMC = ones(length(t_MC), no_runs);
Pureza_MC = ones(length(t_MC), no_runs);
Pureza_isoMC = ones(length(t_MC), no_runs);
Rend_MC = ones(length(t_MC), no_runs);
frac_store = zeros(no_runs, length(frac));
Tstoredin = zeros(no_runs, length(t_MC), N + 1);

tic

for j = 1:no_runs
    % Monte Carlo Loop

    [frac] = getRandFrac();
    [ReactorOut] = getVolumeReactor(frac, TdegC0, spec);
    frac_store(j, :) = frac;

    for i = 1:length(t_MC)
        treaction = t_MC(i); % h

        % Isotérmico
        fun = @(T) OptimizeReactorIso(treaction, frac, spec, T);
        [TstarMC(i, j), obj_isoMC(i, j)] = fminbnd(fun, Tmin, Tmax);
        metricas_TstarMC = OptimizeReactorIso2(treaction, frac, spec, TstarMC(i, j));
        Pureza_isoMC(i, j) = metricas_TstarMC.Pureza;

        % Dinâmico
        [metricas_opti, solT, solC] = RendDyn(N, treaction, TdegC0, 'rk', Tmin, Tmax, ReactorOut.C0_sorted);
        Tstoredin(j, i, :) = solT;
        Pureza_MC(i, j) = metricas_opti.Pureza;
        Rend_MC(i, j) = metricas_opti.PercentagemHidrogenada;
    end

    j
end

time_MonteCarlo = toc;

% Dinâmica
OBJ_MC_mean = 100 * mean(Rend_MC(:, 1:60), 2);
OBJ_MC_std = 100 * std(Rend_MC(:, 1:60), [], 2);

Pureza_MC_mean = 100 * mean(Pureza_MC(:, 1:60), 2);
Pureza_MC_std = 100 * std(Pureza_MC(:, 1:60), [], 2);

% Iso
OBJ_MC_isomean = -100 * mean(obj_isoMC(:, 1:60), 2);
OBJ_MC_isostd = -100 * std(obj_isoMC(:, 1:60), [], 2);

Pureza_MC_isomean = 100 * mean(Pureza_isoMC(:, 1:60), 2);
Pureza_MC_isostd = 100 * std(Pureza_isoMC(:, 1:60), [], 2);

%%
close all

figure(Color = "w")
plot(t_MC, OBJ_MC_mean, LineWidth = 2, DisplayName = "\psi", LineStyle = "-"); hold on
plot(t_MC, OBJ_MC_mean - OBJ_MC_std, LineWidth = 1, DisplayName = "\psi - \sigma", LineStyle = ":", Color = [0.11, 0.62, 0.47]); hold on
plot(t_MC, OBJ_MC_mean + OBJ_MC_std, LineWidth = 1, DisplayName = "\psi + \sigma", LineStyle = ":", Color = [0.11, 0.62, 0.47]); hold on

plot(t_MC, OBJ_MC_isomean, LineWidth = 2, DisplayName = "ISO \psi", LineStyle = "-"); hold on
plot(t_MC, OBJ_MC_isomean - OBJ_MC_isostd, LineWidth = 1, DisplayName = "\psi - \sigma", LineStyle = ":", Color = [0.91, 0.16, 0.54]); hold on
plot(t_MC, OBJ_MC_isomean + OBJ_MC_isostd, LineWidth = 1, DisplayName = "\psi + \sigma", LineStyle = ":", Color = [0.91, 0.16, 0.54]); hold on

plot(t_MC, Pureza_MC_mean, LineWidth = 2, DisplayName = "\eta (%)", LineStyle = "--"); hold on
plot(t_MC, Pureza_MC_mean - Pureza_MC_std, LineWidth = 1, DisplayName = "\eta (%) - \sigma", LineStyle = ":", Marker = "none", Color = [0.65, 0.46, 0.11]); hold on
plot(t_MC, Pureza_MC_mean + Pureza_MC_std, LineWidth = 1, DisplayName = "\eta (%) + \sigma", LineStyle = ":", Marker = 'none', Color = [0.65, 0.46, 0.11]); hold on

plot(t_MC, Pureza_MC_isomean, LineWidth = 2, DisplayName = "ISO \eta (%)", LineStyle = "-.", Color = [0.46, 0.44, 0.70]); hold on
plot(t_MC, Pureza_MC_isomean - Pureza_MC_isostd, LineWidth = 1, DisplayName = "\eta (%) - \sigma", LineStyle = ":", Marker = "none", Color = [0.46, 0.44, 0.70]); hold on
plot(t_MC, Pureza_MC_isomean + Pureza_MC_isostd, LineWidth = 1, DisplayName = "\eta (%) + \sigma", LineStyle = ":", Marker = 'none', Color = [0.46, 0.44, 0.70]); hold on

yline(92, DisplayName = "Referência 1", LineStyle = "--")
yline(90, DisplayName = "Referência 2", LineStyle = ":")

legend(Box = "off")
ylabel("\psi, \eta (%)")
xlabel("t [h]")
set(gca, "FontSize", 12)
