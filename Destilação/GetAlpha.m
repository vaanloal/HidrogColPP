function [alpha_mean, alpha] = GetAlpha(Ant, P, visualize)

    P_vap_Col = @(TdegC) 10 .^ (Ant.Col(1) - Ant.Col(2) ./ (TdegC + Ant.Col(3)));
    P_vap_Xyl = @(TdegC) 10 .^ (Ant.Xil(1) - Ant.Xil(2) ./ (TdegC + Ant.Xil(3)));
    T_vap_Col = @(P) (Ant.Col(2) + log10(P) * Ant.Col(3) - (Ant.Col(1) * Ant.Col(3))) / (Ant.Col(1) - log10(P));
    T_vap_Xyl = @(P) (Ant.Xil(2) + log10(P) * Ant.Xil(3) - (Ant.Xil(1) * Ant.Xil(3))) / (Ant.Xil(1) - log10(P));
    alpha_LK_HK = @(T) P_vap_Xyl(T) ./ P_vap_Col(T);

    x = linspace(0, 1, 100); % fração do LK
    T_mix = x * T_vap_Xyl(P) + (1 - x) * T_vap_Col(P);
    y_xil = x .* P_vap_Xyl(T_mix) / P;
    alpha = alpha_LK_HK(T_mix);
    alpha_mean = mean(alpha);

    if visualize == 1
        hfig = figure(Color = "w");
        plot(T_mix, log10(alpha), 'LineWidth', 1.5)
        xlabel("$T_{\rm mix}\; \rm[^oC]$", Interpreter = "latex")
        ylabel("$\log_{10}(\alpha_{\rm LK/HK})$", Interpreter = "latex")
        ax = gca;
        SetFigureDefaults2(hfig, 0.5, 15, 12, ax)
        axis tight
    else
    end

end
