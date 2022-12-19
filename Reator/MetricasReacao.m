function [metricas] = MetricasReacao(C)
    % AA(1)   DIA(2) TEA(3) PAA(4)  DEA(5) NEA(6)  PIA(7) DPIAs(8)
    OBJ_SIM1 = (C(:, 8) + C(:, 2) + C(:, 3) - (C(1, 8) + C(1, 2) + C(1, 3))) ./ ...
    (C(1, 1) + C(1, 4) + C(1, 6) + C(1, 7));

    xAA = (C(1, 1) - C(end, 1)) / C(1, 1);
    xPIAs = (C(1, 7) - C(end, 7)) / C(1, 7);
    xNEA = (C(1, 6) - C(end, 6)) / C(1, 6);
    xPAA = (C(1, 4) - C(end, 4)) / C(1, 4);
    S_DIA_DEA = C(end, 2) / C(end, 5);
    S_TEA_DEA = C(end, 3) / C(end, 5);
    % S_inst_1 = (C1(:, 3) + C1(:, 2))./C1(:, 5);

    Pureza = (C(end, 1) + C(end, 2) + C(end, 3) + C(end, 4) ...
    + C(end, 6) + C(end, 7) + C(end, 8)) / (C(end, 1) ...
        + C(end, 2) + C(end, 3) + C(end, 4) ...
        + C(end, 6) + C(end, 7) + C(end, 8) + C(end, 5));

    metricas.PercentagemHidrogenada = OBJ_SIM1(end);
    metricas.Pureza = Pureza;
    metricas.xAA = xAA;
    metricas.xPIAS = xPIAs;
    metricas.xNEA = xNEA;
    metricas.xPAA = xPAA;
    metricas.S_DIA_TEA = S_DIA_DEA;
    metricas.S_TEA = S_TEA_DEA;

end
