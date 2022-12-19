function [out] = getVolumeReactor(frac, T, spec)
    % AA(1)   DIA(2) TEA(3) PAA(4)  DEA(5) NEA(6)  PIA(7) DPIAs(8)
    Mw = [302 302 302 304 300 306 302 304];

    HeadSpace = spec.HeadSpace;
    TC = spec.TC; % h
    TLplus = spec.TLplus; % h
    ProducaoAnual = spec.ProducaoAnual; % ton/ano
    SemanasAno = spec.SemanasAno;
    HorasDia = spec.HorasDia; % h
    DiasSemana = spec.DiasSemana;
    HD_racio = spec.HD_racio;

    HorasSemana = HorasDia * DiasSemana;
    HorasAno = HorasSemana * SemanasAno;

    N = (HorasAno - (TC + TLplus)) ./ (TC) - 1;

    ProducaoAnualMolarCH = ProducaoAnual * 1000 * 1000/304.5;
    ReacoesAno = fix(N);
    ProducaoBatch = ProducaoAnualMolarCH / ReacoesAno; % mol/batch

    DIA = ProducaoBatch ./ sum(frac ./ frac(2));
    NEA = (DIA * frac(6)) ./ frac(2);
    PAA = (DIA * frac(4)) ./ frac(2);
    AA = (DIA * frac(1)) ./ frac(2);
    DEA = (DIA * frac(5)) ./ frac(2);
    TEA = (DIA * frac(3)) ./ frac(2);
    PIAs = (DIA * frac(7)) ./ frac(2);
    DPIAs = (DIA * frac(8)) ./ frac(2);
    Xileno_mol = sum([NEA PAA AA DIA DEA TEA PIAs DPIAs] .* Mw) / 106.16; % gramas de colofónia = gramas de xileno
    Xileno_massa = sum([NEA PAA AA DIA DEA TEA PIAs DPIAs] .* Mw) / 1000; % kg

    % Calcular o Volume do Reator. Através das massas volúmicas, da massa de catalisador e do espaço alocado para a serpentina.

    MassaCol = sum([NEA PAA AA DIA DEA TEA PIAs DPIAs] .* Mw) / 1000; % kg
    VolumeCol = MassaCol / rho_col(T); % m3
    VolumeXil = Xileno_massa / rho_xyl(T); % m3

    rhoNiquel = 8.908; % g/cm3
    rhoCarvao = 2.1; % g/cm3
    MassaMetal = sum([NEA PAA AA DIA DEA TEA PIAs DPIAs] .* Mw) * 0.05 * 0.05;
    VolumeMetal = MassaMetal / rhoNiquel; % cm3

    MassaCarvao = sum([NEA PAA AA DIA DEA TEA PIAs DPIAs] .* Mw) * 0.05 * 0.95;
    VolumeCarvao = MassaCarvao / rhoCarvao; % cm3

    VolumeCatalisador = (VolumeCarvao + VolumeMetal) / (10 ^ 6); % m3

    VolumeLiquid = VolumeCatalisador + VolumeXil + VolumeCol; % m3
    VReactor = (HeadSpace + 1) * VolumeLiquid;

    DReactor = (4 * VReactor / (HD_racio * pi)) .^ (1/3); % m
    HReactor = VReactor / (pi * (DReactor / 2) ^ 2); % m
    AreaCamisa = pi * DReactor / 2 * (DReactor / 2 + 2 * HReactor); % m2

    C0 = [NEA PAA AA DIA DEA TEA PIAs DPIAs] ./ (VolumeLiquid * 1000); % mol/L
    C0_sorted = C0([3 4 6 2 5 1 7 8]);

    out.ProducaoBatch = ProducaoBatch;
    out.V_Liquid = VolumeLiquid;
    out.V_Reactor = VReactor;
    out.D_Reactor = DReactor;
    out.H_Reactor = HReactor;
    out.Area_Camisa = AreaCamisa;
    out.C0 = C0;
    out.C0_sorted = C0_sorted;
    out.MassaMix = MassaCol + Xileno_massa;

end
