function [metricas_opti, res] = RendDynTime(frac, TdegC0, spec, N, Tmin, Tmax, tr_min, tr_max, IC, RendimentoGOAL)

    opti = casadi.Opti(); % Optimization problem

    [ReactorOut] = getVolumeReactor(frac, TdegC0, spec);

    x0 = ReactorOut.C0_sorted;

    % Modelo
    state = casadi.MX.sym('states', 8); % 8 states
    controls = casadi.MX.sym('controls', 1); % 1 control
    xdot = rhs_optimize_time(state, controls);
    f = casadi.Function('f', {state, controls}, {xdot});

    X = opti.variable(8, N + 1); % state trajectory
    C_AA = X(1, :);
    C_DIA = X(2, :);
    C_TEA = X(3, :);
    C_PAA = X(4, :);
    C_DEA = X(5, :);
    C_NEA = X(6, :);
    C_PIAs = X(7, :);
    C_DPIAs = X(8, :);

    U = opti.variable(1, N); % control trajectory (temperature)
    Time = opti.variable(); % final time

    opti.minimize(Time);

    dt = Time / N; % length of a control interval

    for k = 1:N % loop over control intervals
        % Runge-Kutta 4 integration
        k1 = f(X(:, k), U(:, k));
        k2 = f(X(:, k) + dt / 2 * k1, U(:, k));
        k3 = f(X(:, k) + dt / 2 * k2, U(:, k));
        k4 = f(X(:, k) + dt * k3, U(:, k));
        x_next = X(:, k) + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4);
        opti.subject_to(X(:, k + 1) == x_next); % close the gaps
    end

    opti.subject_to(Tmin + 273.15 <= U <= Tmax + 273.15);

    opti.subject_to(C_AA(1) == x0(1));
    opti.subject_to(C_DIA(1) == x0(2));
    opti.subject_to(C_TEA(1) == x0(3));
    opti.subject_to(C_PAA(1) == x0(4));
    opti.subject_to(C_DEA(1) == x0(5));
    opti.subject_to(C_NEA(1) == x0(6));
    opti.subject_to(C_PIAs(1) == x0(7));
    opti.subject_to(C_DPIAs(1) == x0(8));

    opti.subject_to(tr_min * 60 <= Time <= tr_max * 60); % Time must be positive

    opti.subject_to(C_AA(1, 2:end) <= 2);
    opti.subject_to(C_DIA(1, 2:end) <= 2);
    opti.subject_to(C_TEA(1, 2:end) <= 2);
    opti.subject_to(C_PAA(1, 2:end) <= 2);
    opti.subject_to(C_DEA(1, 2:end) <= 2);
    opti.subject_to(C_NEA(1, 2:end) <= 2);
    opti.subject_to(C_PIAs(1, 2:end) <= 2);
    opti.subject_to(C_DPIAs(1, 2:end) <= 2);

    opti.set_initial(C_AA(1, 2:end), IC(:, 1));
    opti.set_initial(C_DIA(1, 2:end), IC(:, 2));
    opti.set_initial(C_TEA(1, 2:end), IC(:, 3));
    opti.set_initial(C_PAA(1, 2:end), IC(:, 4));
    opti.set_initial(C_DEA(1, 2:end), IC(:, 5));
    opti.set_initial(C_NEA(1, 2:end), IC(:, 6));
    opti.set_initial(C_PIAs(1, 2:end), IC(:, 7));
    opti.set_initial(C_DPIAs(1, 2:end), IC(:, 8));
    opti.set_initial(U, TdegC0 + 273.15);
    opti.set_initial(Time, 60); % min

    % Pureza = (C_AA(end) + C_DIA(end) + C_TEA(end) + C_PAA(end) + ...
    %     C_NEA(end) + C_PIAs(end) + C_DPIAs(end)) / (C_AA(end) + ...
    %     C_DIA(end) + C_TEA(end) + C_PAA(end) + C_DEA(end) + ...
    %     C_NEA(end) + C_PIAs(end) + C_DPIAs(end));

    Rend = (C_DIA(end) + C_TEA(end) + C_DPIAs(end) - (x0(8) + x0(2) + x0(3))) ./ (x0(1) + x0(4) + x0(6) + x0(7));

    opti.subject_to(Rend == RendimentoGOAL);

    opti.solver('ipopt', struct('print_time', false), struct('print_level', 0)); % set numerical backend
    % opti.callback(@(i) plot(opti.debug.value(X')))
    sol = opti.solve(); % actual solve

    res.states = opti.value(X)';
    res.controls = opti.value(U)' - 273.15; %oC
    res.dt = opti.value(dt) / 60; % h
    res.tf = opti.value(Time) / 60; % h
    res.tgrid = 0:res.dt:res.tf;
    res.FinalTime = opti.value(Time) / 60; % h

    metricas_opti = MetricasReacao(res.states);

end
