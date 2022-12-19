clc; clear all; close all; format long

% Adicionar o CasADi ao path do Matlab
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');

load("500_MC_min_time.mat")
clearvars -except res metricas_opti

% Exemplo 1
i = 9;

for i = 1:500

    T_profile = [res(i).controls; res(i).controls(end)]; % perfil do modelo dinâmico

    T_slope = (res(i).controls(end) - res(i).controls(1)) ./ ...
        ((res(i).tgrid(end) - (res(i).tgrid(1))) * 60); % [oC/min]

    T_linear_fun = @(tmin) T_slope .* tmin + T_profile(1); % função para calcular o perfil
    tgrid = res(i).tgrid'; % intervalos de tempo do modelo dinâmico

    C_time_opti = res(i).states; % armazenar perfis das variáveos de estado otimizadas
    C0 = res(i).states(1, :); % condições iniciais para o modelo novo

    % INPUT Function
    time_min = tgrid * 60; % min
    input_T = T_linear_fun(time_min); % perfil linear de Temperatura [oC]
    % input_T(end) = nan; % último valor de Temperatura não afeta os estados

    % Simular o novo perfil de temperaturas
    odeFun = @(t, x) rhs_linearT(t, x, input_T, time_min);
    [t_slopeT, C_slopeT] = ode15s(odeFun, time_min, C0');
    [metricas_slopeT] = MetricasReacao(C_slopeT);

    PurezaLinear(i) = metricas_slopeT.Pureza;
    RendimLinear(i) = metricas_slopeT.PercentagemHidrogenada;

    PurezaOpti(i) = metricas_opti(i).Pureza;
    RendimOpti(i) = metricas_opti(i).PercentagemHidrogenada;

    profileT_linear(:, i) = input_T;
    profileT_opti(:, i) = T_profile;

end

%% PLot Temperature Profiles
close all

figure(Color = "W", Position = [680, 663, 742, 315])
t = tiledlayout(1, 2);
t.TileSpacing = 'compact';
t.Padding = 'compact';
nexttile

for i = 1:500
    plot(res(i).tgrid, profileT_opti(:, i)); hold on
end

xlabel("$t\; \rm[h]$", Interpreter = "latex")
ylabel("$T\; \rm[^oC]$", Interpreter = "latex")
text(0, 200, '$a)$', Interpreter = 'latex', FontSize = 12)
set(gca, FontSize = 12) %, PlotBoxAspectRatio = [1,1,1])
axis padded

nexttile

for i = 1:500
    plot(res(i).tgrid, profileT_linear(:, i)); hold on
end

xlabel("$t\; \rm[h]$", Interpreter = "latex")
ylabel("$T\; \rm[^oC]$", Interpreter = "latex")
text(0, 200, '$b)$', Interpreter = 'latex', FontSize = 12)
set(gca, FontSize = 12) %, PlotBoxAspectRatio = [1,1,1])
axis padded

% exportgraphics(gcf, 'temp_profiles.pdf', Resolution=600)
%% Gráfico da Pureza e do Rendimento
close all
figure(Color = "W", Position = [680, 663, 742, 315])
t = tiledlayout(1, 2);
t.TileSpacing = 'compact';
t.Padding = 'compact';

nexttile
plot(RendimOpti * 100, LineWidth = 1.5, LineStyle = "--", DisplayName = 'Dinâmico'); hold on
plot(RendimLinear * 100, LineWidth = 1, DisplayName = 'Rampa'); hold on
set(gca, "FontSize", 12)
xlabel('Simulação', FontName = 'Georgia')
ylabel("$\psi\; [\%]$", Interpreter = "latex")
ylim([93, 95])
legend(Box = "off", FontSize = 12, FontName = 'Georgia')

nexttile
A = {PurezaOpti' * 100, PurezaLinear' * 100};
nhist(A, 'legend', {'Dinâmico', 'Rampa'});

% t1 = nhist([PurezaOpti'*100, PurezaLinear'*100],'pdf','legend',{'Dinâmico','Rampa'},'box'); hold on
% t1 = nhist(PurezaOpti'*100,'legend','Dinâmico', 'linewidth',1,'box','color',[0.11,0.62,0.47]); hold on
% t2 = nhist(PurezaLinear'*100,'legend','Rampa', 'linewidth',1,'box','color',[0.85,0.37,0.01]); hold on
% face_alpha = 0.2;
% edge_alpha = 0.2;
%
% edge_color = [0 0 0];%[0.96 0.96 0.96];
% face_color = [0.11,0.62,0.47];
% linwidth = 0.3;
% methods = {'fd'; 'scott'; 'sturges'};
% idx = 3;
% n = calcnbins(PurezaOpti*100, methods{idx});
%
% h1 = histfit(PurezaOpti*100,n,'kernel', FaceColor= face_color);
% set(h1(1),'FaceAlpha',face_alpha);
% set(h1(1),'FaceColor',face_color)
% set(h1(1),'EdgeColor',edge_color);
% set(h1(1),'LineWidth',linwidth);
% xlabel('$t\;\rm [h]$', Interpreter='latex')
% ylabel("$\rm Frequ\hat{e}ncia$", Interpreter="latex")
% set(gca, "FontSize",12)
% hold on
%
%
% edge_color = [0 0 0];%[0.96 0.96 0.96];
% face_color = [0.85,0.37,0.01];
% n = calcnbins(PurezaLinear*100, methods{idx});
% h1 = histfit(PurezaLinear*100,n,'kernel', FaceColor= face_color);
% set(h1(1),'FaceAlpha',face_alpha);
% set(h1(1),'FaceColor',face_color)
% set(h1(1),'EdgeColor',edge_color);
% set(h1(1),'LineWidth',linwidth);
% xlabel('$\chi\;\rm [\%]$', Interpreter='latex')
% ylabel("$\rm Frequ\hat{e}ncia$", Interpreter="latex")
% set(gca, "FontSize",12)
%
% % figure(Color="W")
% % plot(tgrid, T_profile, LineWidth=2); hold on
% % plot(tgrid, input_T, LineWidth=2); hold on
% % plot(t_slopeT/60, input_T, LineWidth=2); hold on
% %
% % xlabel("t [h]")
% % ylabel("T [^oC]")

%% Simulação da Rampa de T com quenching
%%%%%%%% MAIN INPUTS %%%%%%%%
treaction = 2.7; % h
toutros = 2; % h
TdegC0 = 200; % oC
Tfiltro_halt_degC = 150; % oC
Tcoolant_degC = mean([20, 50]); % oC
U_cooling = 400; % 500 - 1200 W/m2/K
tempo_transicao = 5; % min

%%%%%% END MAIN INPUTS %%%%%%

TK0 = TdegC0 + 273.15; % K
% frac = [0.487772819 0.046512298 0.028606858 0.164908937 0.035763384 0.138086399 0.098349306 0];
frac = [0.309422684289359	0.0458880722748179	0.0625631193684445	0.286655191060494	0.0583190052412689	0.145108410333857	0.0895904633869377	0.00245305404482194];

spec.HeadSpace = 0.2; %
spec.TC = treaction + toutros; % h
spec.TLplus = 13; % h
spec.ProducaoAnual = 4785; % ton/ano
spec.SemanasAno = 48;
spec.HorasDia = 24; % h
spec.DiasSemana = 7;
spec.HD_racio = 1.3; % Razão Altura/Diâmetro do reator

[ReactorOut] = getVolumeReactor(frac, TdegC0, spec);

% Simulação Inicial
t_span = linspace(0, treaction * 60, 300);
T_slope = (mean(profileT_linear(end, :)) - mean(profileT_linear(1, :))) ./ ...
    (t_span(end) - t_span(1)); % [oC/min]

T_linear_fun = @(tmin) T_slope .* tmin + mean(profileT_linear(1, :)); % função para calcular o perfil

% INPUT Function
input_T = T_linear_fun(t_span); % perfil linear de Temperatura [oC]
% input_T(end) = nan; % último valor de Temperatura não afeta os estados

% Simular o novo perfil de temperaturas
odeFun = @(t, x) rhs_linearT(t, x, input_T, t_span);
[t1, C1] = ode15s(odeFun, t_span, ReactorOut.C0_sorted');
[metricas1] = MetricasReacao(C1);

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
T_all = [input_T' + 273.15; C2(2:end, end)] - 273.15;

% Métricas da Simulação II
[metricas_simII] = MetricasReacao(C_all);

%% PLOT KINETICS
nomes = ["$C_{AA}$", "$C_{DIA}$", "$C_{TEA}$", "$C_{PAA}$", "$C_{DEA}$", "$C_{NEA}$", "$C_{PIAs}$", "$C_{DPIAs}$"];
all_marks = {'o', '+', '*', 'h', 'x', 's', 'd', '^', 'v', '>', '<', 'p', '.'};
hfig = figure("Color", "w", Position = [150, 241.5, 654.75, 281.25]);
t = tiledlayout(1, 2)
t.TileSpacing = 'compact';
t.Padding = 'compact';
lin_width = 1;
nexttile
hold on
npoints = length(C_all);
indicis_markers = [10:30:npoints];
plot(t_all / 60, C_all(:, 1), LineWidth = lin_width, DisplayName = "$C_{AA}$", Marker = all_marks{1}, MarkerIndices = indicis_markers)
plot(t_all / 60, C_all(:, 2), LineWidth = lin_width, DisplayName = "$C_{DIA}$", Marker = all_marks{2}, MarkerIndices = indicis_markers)
plot(t_all / 60, C_all(:, 3), LineWidth = lin_width, DisplayName = "$C_{TEA}$", Marker = all_marks{3}, MarkerIndices = indicis_markers)
plot(t_all / 60, C_all(:, 4), LineWidth = lin_width, DisplayName = "$C_{PAA}$", Marker = all_marks{4}, MarkerIndices = indicis_markers)
plot(t_all / 60, C_all(:, 5), LineWidth = lin_width, DisplayName = "$C_{DEA}$", Marker = all_marks{5}, MarkerIndices = indicis_markers)
plot(t_all / 60, C_all(:, 6), LineWidth = lin_width, DisplayName = "$C_{NEA}$", Marker = all_marks{6}, MarkerIndices = indicis_markers)
plot(t_all / 60, C_all(:, 7), LineWidth = lin_width, DisplayName = "$C_{PIAs}$", Marker = all_marks{7}, MarkerIndices = indicis_markers)
plot(t_all / 60, C_all(:, 8), LineWidth = lin_width, DisplayName = "$C_{DPIAs}$", Marker = all_marks{8}, MarkerIndices = indicis_markers)
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

% vec_excel = [xAA; xPIAs; xNEA;xPAA;S_DIA_DEA;S_TEA_DEA];

vec_excel = [metricas_simII.xAA; metricas_simII.xPIAS; metricas_simII.xNEA; metricas_simII.xPAA; metricas_simII.S_DIA_TEA; metricas_simII.S_TEA]
