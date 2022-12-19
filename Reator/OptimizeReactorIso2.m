function [metricas_simI] = OptimizeReactorIso2(treaction, frac, spec, TdegC0)

    [ReactorOut] = getVolumeReactor(frac, TdegC0, spec);

    % Simulação Inicial
    odeFun1 = @(t, x) rhscorr(t, x, TdegC0 + 273.15);
    [~, C] = ode15s(odeFun1, [0, treaction * 60], ReactorOut.C0_sorted');

    % Métricas da Simulação I
    [metricas_simI] = MetricasReacao(C);
end
