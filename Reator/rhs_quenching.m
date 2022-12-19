function xdot = rhs_quenching(t, x, par)

    Tcoolant = par(1);
    m_mix = par(2);
    cp_mix = par(3);
    U = par(4);
    A = par(5);
    tempo_transicao = par(6);
    T_filtro = par(7);

    C_AA = x(1);
    C_DIA = x(2);
    C_TEA = x(3);
    C_PAA = x(4);
    C_DEA = x(5);
    C_NEA = x(6);
    C_PIAs = x(7);
    C_DPIAs = x(8);
    T = x(9); % K

    k1 = 4.54e6 * exp(-66.60 * 1000 / (8.3145 * T));
    k2 = 1.51e5 * exp(-58.06 * 1000 / (8.3145 * T));
    k3 = 1.19e7 * exp(-73.02 * 1000 / (8.3145 * T));
    k4 = 2.35e8 * exp(-86.67 * 1000 / (8.3145 * T));
    k5 = 1.43e10 * exp(-107.42 * 1000 / (8.3145 * T));
    k6 = 2e4 * exp(-37.60 * 1000 / (8.3145 * T));

    r1 = k1 * C_NEA;
    r2 = k2 * C_PAA;
    r3 = k3 * C_AA;
    r4 = k4 * C_AA;
    r5 = k5 * C_AA;
    r6 = k6 * C_PIAs;

    dC_AA =+ r1 + r2 - r3 - r4 - r5;
    dC_DIA =+ r3;
    dC_TEA =+ r4;
    dC_PAA =- r2;
    dC_DEA =+ r5;
    dC_NEA =- r1;
    dC_PIAs =- r6;
    dC_DPIAs =+ r6;

    if t >= tempo_transicao % Nos primeiros 5 min ainda ocorre reação

        dC_AA = 0;
        dC_DIA = 0;
        dC_TEA = 0;
        dC_PAA = 0;
        dC_DEA = 0;
        dC_NEA = 0;
        dC_PIAs = 0;
        dC_DPIAs = 0;

    end

    dTdt = 60 * U * A * (Tcoolant - T) / (m_mix * cp_mix);

    xdot = [dC_AA; dC_DIA; dC_TEA; dC_PAA; dC_DEA; dC_NEA; dC_PIAs; dC_DPIAs; dTdt];
end
