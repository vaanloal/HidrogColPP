function [metricas_opti, solT, solC] = RendDyn(N, treaction, TdegC0, choose_solver, Tmin, Tmax, C0)

    % N = 150;
    % choose_solver = 'rk';

    time = treaction * 60; % tempo de simulação [min]
    tgrid = linspace(0, time, N + 1);

    % Declare model variables
    C_AA = casadi.MX.sym('C_AA');
    C_DIA = casadi.MX.sym('C_DIA');
    C_TEA = casadi.MX.sym('C_TEA');
    C_PAA = casadi.MX.sym('C_PAA');
    C_DEA = casadi.MX.sym('C_DEA');
    C_NEA = casadi.MX.sym('C_NEA');
    C_PIAs = casadi.MX.sym('C_PIAs');
    C_DPIAs = casadi.MX.sym('C_DPIAs');
    x = [C_AA; C_DIA; C_TEA; C_PAA; C_DEA; C_NEA; C_PIAs; C_DPIAs];

    % Controls
    % Temperature (K)
    T = casadi.MX.sym('T');
    T_lb = Tmin + 273.15; % [K]
    T_ub = Tmax + 273.15; % [K]
    T0 = TdegC0 + 273.15; % [K]
    u = T;
    u0 = T0;

    % Model equations
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

    dC_AA = r1 + r2 - r3 - r4 - r5;
    dC_DIA = r3;
    dC_TEA = r4;
    dC_PAA =- r2;
    dC_DEA = r5;
    dC_NEA =- r1;
    dC_PIAs =- r6;
    dC_DPIAs = r6;

    xdot = [dC_AA; dC_DIA; dC_TEA; dC_PAA; dC_DEA; dC_NEA; dC_PIAs; dC_DPIAs];
    F = casadi.Function('F', {x, u}, {xdot}, {'x', 'u'}, {'xdot'});

    % Integrador: discretizar o sistema
    intg_options = struct;

    if choose_solver == 'rk'
        intg_options.tf = time / N;
        intg_options.simplify = true;
        intg_options.number_of_finite_elements = 4;
    else
        intg_options.tf = time / N;
    end

    % DAE problem structure
    dae = struct;
    dae.x = x; % estados do sistema
    dae.p = u; % parâmetros
    dae.ode = F(x, u);

    intg = casadi.integrator('intg', choose_solver, dae, intg_options);
    res = intg('x0', x, 'p', u); x_next = res.xf;

    % Simplificar a API para: (x,u) -> (x_next)
    fun = casadi.Function('fun', {x, u}, {x_next}, {'x', 'u'}, {'x_next'});
    sim = fun.mapaccum(N); % chamar [fun] recursivamente

    x0 = C0';

    resSim = full([x0, sim(x0, u0)])';
    RendSim1 = (resSim(:, 8) + resSim(:, 2) + resSim(:, 3) - (resSim(1, 8) + resSim(1, 2) + resSim(1, 3))) ./ ...
        (resSim(1, 1) + resSim(1, 4) + resSim(1, 6) + resSim(1, 7));
    RendSim1(end)

    %% Otimização dinâmica
    opti = casadi.Opti();

    % determinar as variáveis de estado (quais e quantas)
    % escalar variáveis
    x = opti.variable(8, N + 1);
    C_AA = x(1, :);
    C_DIA = x(2, :);
    C_TEA = x(3, :);
    C_PAA = x(4, :);
    C_DEA = x(5, :);
    C_NEA = x(6, :);
    C_PIAs = x(7, :);
    C_DPIAs = x(8, :);

    % controlos
    T = opti.variable(1, N);
    u = T;

    p = opti.parameter(size(x, 1), 1); %  permanecem fixos durante a simulação

    % função objetivo
    obj = (C_DIA(end) + C_TEA(end) + C_DPIAs(end) - (x0(8) + x0(2) + x0(3))) ./ (x0(1) + x0(4) + x0(6) + x0(7));

    opti.minimize(-obj)

    % Inicializações e restrições do multiple shooting
    for k = 1:N
        opti.subject_to(x(:, k + 1) == fun(x(:, k), u(:, k)));

        % Restrição variáveis de controlo
        opti.subject_to(T_lb <= u(1, k) <= T_ub);
        opti.set_initial(u(1, k), T0);

        % inicializar as variáveis
        for j = 1:size(x, 1)
            opti.set_initial(x(j, k + 1), resSim(k + 1, j));
        end

    end

    % especificar as condições iniciais
    opti.set_value(p, x0);
    % aplicar as condições iniciais
    opti.subject_to(x(:, 1) == p);

    % Limitar a quantidade de ácidos desidrogenados

    % Escolher o solver
    opti.solver('ipopt', struct('print_time', false), struct('print_level', 0));

    % disp("Waiting...")
    % tic
    optisol = opti.solve();
    % toc
    % disp("Solved!")

    solC = opti.value(x)';
    solT = opti.value([u nan])' - 273.15;

    [metricas_opti] = MetricasReacao(solC);

end
