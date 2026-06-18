%% ===========================================================================
%  Phase 1: CSTR Modeling, Linearization, and Sensitivity Analysis
%  Non-isothermal Continuous Stirred Tank Reactor (CSTR)
%  Based on the dimensionless model from Perez & Albertos (2004)
% ===========================================================================
%  States    : x = [x2; x3]  (dimensionless concentration and temperature)
%  Input     : u = x5        (dimensionless cooling water flow rate)
%  Output    : y = x3        (dimensionless reactor temperature)
% ===========================================================================
%  This script performs:
%    1. Parameter definition and dimensionless constants
%    2. Nonlinear state-space model (symbolic)
%    3. Equilibrium point calculation at nominal input
%    4. Linearization around the equilibrium (Jacobians)
%    5. Stability analysis (poles, damping, step response)
%    6. Sensitivity analysis (±10% change in heat transfer coefficient U)
%    7. Canonical forms (controller, observer, controllability, observability)
%    8. Jordan form
%    9. Comprehensive controllability & observability analysis (5 methods)
%   10. Model validation: nonlinear vs. linear comparison (step & sinusoidal)
%   11. RMSE analysis vs. input amplitude and frequency
% ===========================================================================
clear; clc; close all;

%% ===========================================================================
%  SECTION 1: Physical Constants and Dimensionless Parameters
% ===========================================================================
fprintf('========== 1. PARAMETERS ==========\n');

% --- Physical constants (Table 1 from reference) ---
V_s    = 1.36;      % m^3, steady state reactor volume
F_os   = 1.13;      % m^3/h, nominal inlet volumetric flow rate
C_aor  = 3.92;      % kmol/m^3, initial reactant concentration
C_ao   = 8;         % kmol/m^3, inlet reactant concentration
alpha  = 7.08e10;   % 1/h, pre-exponential factor
H_rxn  = 69815;     % kJ/kmol, enthalpy of reaction
E_a    = 69815;     % kJ/kmol, activation energy (same as H_rxn)
R_gas  = 8.314;     % kJ/(kmol·K), perfect gas constant
rho    = 800;       % kg/m^3, density of inlet/outlet streams
Cp     = 3.13;      % kJ/(kg·K), heat capacity of streams
U      = 3065;      % kJ/(h·m^2·K), overall heat transfer coefficient
A_ht   = 23.22;     % m^2, heat transfer area
F_js   = 1.4130;    % m^3/h, steady state cooling water flow rate
V_j    = 0.085;     % m^3, jacket volume
rho_j  = 1000;      % kg/m^3, density of cooling water
Cp_j   = 4.18;      % kJ/(kg·K), heat capacity of cooling water
T_o    = 294.7;     % K, inlet stream temperature
T_jo   = 294.7;     % K, inlet cooling water temperature
T_r    = 309.9;     % K, set point temperature (reference)

% --- Dimensionless parameters (Table 2) ---
c0 = V_s * alpha / F_os;
% c1 = (V_s * alpha * H_rxn * R_gas * C_aor) / (F_os * rho * Cp * E_a)
% Since H_rxn = E_a, they cancel, leaving:
c1 = (V_s * alpha * R_gas * C_aor) / (F_os * rho * Cp);
c2 = U * A_ht / (rho * Cp * F_os);
c3 = V_s * F_js / (F_os * V_j);
c4_star = (rho * Cp * V_s) / (rho_j * Cp_j * V_j);
c4 = c4_star * c2;                     % c4 = c4* * c2

% Dimensionless constants / disturbances
x1       = 1;                          % V/V_s = 1 (constant volume)
x20      = C_ao / C_aor;               % dimensionless inlet concentration
x30      = R_gas * T_o / E_a;          % dimensionless inlet temperature
x40      = R_gas * T_jo / E_a;         % dimensionless inlet cooling water temp.
x60_nom  = 1;                          % F_o / F_os = 1 (nominal inlet flow)

fprintf('c0 = %.6f, c1 = %.6f, c2 = %.6f, c3 = %.6f, c4 = %.6f\n', c0,c1,c2,c3,c4);
fprintf('x20 = %.6f, x30 = %.6f, x40 = %.6f\n\n', x20,x30,x40);

%% ===========================================================================
%  SECTION 2: Nonlinear State-Space Model (symbolic)
% ===========================================================================
fprintf('========== 2. NONLINEAR MODEL ==========\n');

syms x2 x3 u real
% State derivatives (dimensionless time τ)
f1 = (x60_nom/x1)*(x20 - x2) - c0*x2*exp(-1/x3);
f2 = (x60_nom/x1)*(x30 - x3) + c1*x2*exp(-1/x3) - ...
     (c2*c3*u*(x3 - x40)) / (x1*(c3*u + c4));

f = [f1; f2];          % state equations f(x,u)
h = x3;                % output equation y = h(x) = x3

disp('State equations:');
pretty(f1), pretty(f2)

%% ===========================================================================
%  SECTION 3: Equilibrium Point at Nominal Input (u = 1)
% ===========================================================================
fprintf('========== 3. EQUILIBRIUM POINT ==========\n');

u_nom = 1;                         % nominal cooling water flow (dimensionless)
% Initial guess for equilibrium (approximate)
x3_guess = R_gas * T_r / E_a;      % around 0.0369
x2_guess = 0.8;
x0_guess = [x2_guess; x3_guess];

% Solve f(x, u_nom) = 0 numerically
options = optimoptions('fsolve', 'Display', 'off', 'Algorithm', 'trust-region-dogleg');
eqns = @(x) [...
    (x60_nom/x1)*(x20 - x(1)) - c0*x(1)*exp(-1/x(2)); ...
    (x60_nom/x1)*(x30 - x(2)) + c1*x(1)*exp(-1/x(2)) - ...
        (c2*c3*u_nom*(x(2) - x40)) / (x1*(c3*u_nom + c4)) ];
x_ss = fsolve(eqns, x0_guess, options);

fprintf('Equilibrium states at u = 1:\n');
fprintf('  x2_ss = %.6f  (dim. concentration)\n', x_ss(1));
fprintf('  x3_ss = %.6f  (dim. temperature)\n', x_ss(2));
fprintf('  Corresponding physical temperature = %.3f K\n\n', x_ss(2)*E_a/R_gas);

%% ===========================================================================
%  SECTION 4: Linearization – Jacobians and State-Space Matrices
% ===========================================================================
fprintf('========== 4. LINEARIZATION ==========\n');

% Symbolic Jacobians
A_sym = jacobian(f, [x2; x3]);
B_sym = jacobian(f, u);

% Substitute equilibrium values and nominal input
A = double(subs(A_sym, [x2, x3, u], [x_ss(1), x_ss(2), u_nom]));
B = double(subs(B_sym, [x2, x3, u], [x_ss(1), x_ss(2), u_nom]));

% Output matrices: y = x3 => C = [0, 1], D = 0
C = [0, 1];
D = 0;

fprintf('Linearized state-space matrices:\n');
disp('A = '); disp(A);
disp('B = '); disp(B);
disp('C = '); disp(C);
disp('D = '); disp(D);

%% ===========================================================================
%  SECTION 5: Stability Analysis – Poles, Damping, Step Response
% ===========================================================================
fprintf('\n========== 5. STABILITY ANALYSIS ==========\n');

sys_lin = ss(A, B, C, D);
poles = eig(A);
fprintf('Open-loop poles (eigenvalues of A):\n');
disp(poles);

if all(real(poles) < 0)
    fprintf('System is STABLE at the operating point.\n\n');
else
    fprintf('System is UNSTABLE at the operating point.\n\n');
end

% Display natural frequencies and damping
fprintf('Natural frequency and damping for each pole:\n');
damp(sys_lin);

% Pole-zero map
figure;
pzmap(sys_lin);
title('Pole-Zero Map of Linearized CSTR');
grid on; axis equal;

% Open-loop step response
figure;
step(sys_lin);
title('Open-loop Step Response (from u to y = x_3)');
xlabel('Dimensionless time \tau');
ylabel('x_3 (dimensionless temperature)');
grid on;

%% ===========================================================================
%  SECTION 6: Sensitivity Analysis – Variation of Heat Transfer Coefficient U
% ===========================================================================
fprintf('\n========== 6. SENSITIVITY ANALYSIS (U \26110%%) ==========\n');

perturbation = 0.10;                  % 10%
U_values = [U*(1-perturbation), U, U*(1+perturbation)];
colors = {'b', 'k', 'r'};
legend_str = cell(1,3);

% Pre-allocate
pole_data = zeros(2,3);

% Loop over U values
figure; hold on;
for i = 1:3
    U_i = U_values(i);
    % Recalculate parameters depending on U
    c2_i = U_i * A_ht / (rho * Cp * F_os);
    c4_i = c4_star * c2_i;

    % Solve equilibrium for new U
    eqns_i = @(x) [...
        (x60_nom/x1)*(x20 - x(1)) - c0*x(1)*exp(-1/x(2)); ...
        (x60_nom/x1)*(x30 - x(2)) + c1*x(1)*exp(-1/x(2)) - ...
            (c2_i*c3*u_nom*(x(2) - x40)) / (x1*(c3*u_nom + c4_i)) ];
    x_ss_i = fsolve(eqns_i, x0_guess, options);

    % Linearize at new equilibrium
    A_i = double(subs(A_sym, [x2, x3, u, c2, c4], ...
                        [x_ss_i(1), x_ss_i(2), u_nom, c2_i, c4_i]));
    B_i = double(subs(B_sym, [x2, x3, u, c2, c4], ...
                        [x_ss_i(1), x_ss_i(2), u_nom, c2_i, c4_i]));

    poles_i = eig(A_i);
    pole_data(:,i) = poles_i;

    % Plot poles
    plot(real(poles_i), imag(poles_i), 'x', 'MarkerSize', 10, ...
         'LineWidth', 2, 'Color', colors{i});
    legend_str{i} = sprintf('U = %.1f  (%+.0f%%)', U_i, (U_i-U)/U*100);
end
xlabel('Real Axis'); ylabel('Imaginary Axis');
title('Pole Migration due to \26110% Change in U');
legend(legend_str, 'Location', 'best'); grid on; axis equal;

% Print pole positions
fprintf('\nPole positions:\n');
fprintf('U = %.1f: poles = %.4f \261 %.4fj\n', U_values(1), ...
        real(pole_data(1,1)), abs(imag(pole_data(1,1))));
fprintf('U = %.1f: poles = %.4f \261 %.4fj  (nominal)\n', U_values(2), ...
        real(pole_data(1,2)), abs(imag(pole_data(1,2))));
fprintf('U = %.1f: poles = %.4f \261 %.4fj\n', U_values(3), ...
        real(pole_data(1,3)), abs(imag(pole_data(1,3))));

% Step response comparison for different U
figure; hold on;
for i = 1:3
    U_i = U_values(i);
    c2_i = U_i * A_ht / (rho * Cp * F_os);
    c4_i = c4_star * c2_i;
    eqns_i = @(x) [...
        (x60_nom/x1)*(x20 - x(1)) - c0*x(1)*exp(-1/x(2)); ...
        (x60_nom/x1)*(x30 - x(2)) + c1*x(1)*exp(-1/x(2)) - ...
            (c2_i*c3*u_nom*(x(2) - x40)) / (x1*(c3*u_nom + c4_i)) ];
    x_ss_i = fsolve(eqns_i, x0_guess, options);
    A_i = double(subs(A_sym, [x2, x3, u, c2, c4], ...
                        [x_ss_i(1), x_ss_i(2), u_nom, c2_i, c4_i]));
    B_i = double(subs(B_sym, [x2, x3, u, c2, c4], ...
                        [x_ss_i(1), x_ss_i(2), u_nom, c2_i, c4_i]));
    sys_i = ss(A_i, B_i, C, D);
    step(sys_i);
end
title('Open-loop Step Response for Different U Values');
xlabel('\tau'); ylabel('x_3');
legend(legend_str, 'Location', 'best'); grid on;

%% ===========================================================================
%  SECTION 7: Canonical Forms and Jordan Form
% ===========================================================================
fprintf('\n========== 7. CANONICAL FORMS & JORDAN ==========\n');

% Transfer function coefficients: G(s) = (b1 s + b0) / (s^2 + a1 s + a0)
[num, den] = ss2tf(A, B, C, D);
a0 = den(3); a1 = den(2);
b0 = num(3); b1 = num(2);   % Note: num has three elements for SISO, 1st is 0
fprintf('Transfer function: (%g s + %g) / (s^2 + %g s + %g)\n', b1, b0, a1, a0);

% --- Controller Canonical Form (phase-variable) ---
A_ccf = [0, 1; -a0, -a1];
B_ccf = [0; 1];
C_ccf = [b0, b1];
D_ccf = D;
fprintf('\nController Canonical Form:\n');
disp('A = '); disp(A_ccf);
disp('B = '); disp(B_ccf);
disp('C = '); disp(C_ccf);

% --- Observer Canonical Form (dual) ---
A_ocf = A_ccf';
B_ocf = C_ccf';
C_ocf = B_ccf';
D_ocf = D;
fprintf('\nObserver Canonical Form:\n');
disp('A = '); disp(A_ocf);
disp('B = '); disp(B_ocf);
disp('C = '); disp(C_ocf);

% --- Markov parameters for observability/controllability canonical forms ---
M = [1, 0; a1, 1];
g = M \ [b1; b0];          % [g1; g0]
g1 = g(1); g0 = g(2);

% --- Observability Canonical Form (standard) ---
A_obs = [0, 1; -a0, -a1];
B_obs = [g0; g1];
C_obs = [1, 0];
D_obs = D;
fprintf('\nObservability Canonical Form:\n');
disp('A = '); disp(A_obs);
disp('B = '); disp(B_obs);
disp('C = '); disp(C_obs);

% --- Controllability Canonical Form (standard) ---
A_con = [0, -a0; 1, -a1];
B_con = [1; 0];
C_con = [g1, g0];           % order: g1 then g0 as requested
D_con = D;
fprintf('\nControllability Canonical Form:\n');
disp('A = '); disp(A_con);
disp('B = '); disp(B_con);
disp('C = '); disp(C_con);

% --- Jordan Form (using symbolic) ---
A_sym_mat = sym(A);
[V, J_sym] = jordan(A_sym_mat);
J = double(J_sym);
V_mat = double(V);
B_jordan = V_mat \ B;
C_jordan = C * V_mat;
fprintf('\nJordan Form:\n');
disp('J = '); disp(J);
disp('B_j = '); disp(B_jordan);
disp('C_j = '); disp(C_jordan);

%% ===========================================================================
%  SECTION 8: Controllability Analysis (5 Methods)
% ===========================================================================
fprintf('\n========== 8. CONTROLLABILITY ANALYSIS ==========\n');

n = length(A);

%% Method 1 – Controllability Matrix
Co = ctrb(A, B);
fprintf('\n--- Method 1: Controllability Matrix ---\n');
disp('Co = [B, AB] ='); disp(Co);
fprintf('Rank = %d -> %s\n', rank(Co), iif(rank(Co)==n, 'CONTROLLABLE', 'NOT controllable'));

%% Method 2 – Manual Controllability Matrix
AB = A * B;
Co_hand = [B, AB];
fprintf('\n--- Method 2: Manual Controllability Matrix ---\n');
disp('Co_hand = [B, AB] ='); disp(Co_hand);
fprintf('Rank = %d -> %s\n', rank(Co_hand), iif(rank(Co_hand)==n, 'CONTROLLABLE', 'NOT controllable'));

%% Method 3 – Jordan Form Test
fprintf('\n--- Method 3: Jordan Form Test ---\n');
disp('J = '); disp(J);
disp('B_j = '); disp(B_jordan);
fprintf('All rows of B_j corresponding to each Jordan block must be non-zero.\n');
% (visual inspection)

%% Method 4 – Controllability Gramian
fprintf('\n--- Method 4: Controllability Gramian ---\n');
if all(real(eig(A)) < 0)
    Wc = gram(sys_lin, 'c');
    disp('Wc = '); disp(Wc);
    disp('eig(Wc) = '); disp(eig(Wc));
    if all(eig(Wc) > 0)
        fprintf('Positive definite -> CONTROLLABLE\n');
    else
        fprintf('Not positive definite -> NOT controllable\n');
    end
else
    fprintf('System unstable, Gramian not defined.\n');
end

%% Method 5 – Popov-Belevitch-Hautus (PBH) Test
fprintf('\n--- Method 5: PBH Test ---\n');
eigvals = eig(A);
for i = 1:length(eigvals)
    lambda = eigvals(i);
    H = [A - lambda*eye(n), B];
    fprintf('λ = %.4f + %.4fi, rank([A-λI, B]) = %d\n', ...
            real(lambda), imag(lambda), rank(H));
end
fprintf('All ranks must be n (= %d) for controllability.\n', n);

%% ===========================================================================
%  SECTION 9: Observability Analysis (5 Methods)
% ===========================================================================
fprintf('\n========== 9. OBSERVABILITY ANALYSIS ==========\n');

%% Method 1 – Observability Matrix
Ob = obsv(A, C);
fprintf('\n--- Method 1: Observability Matrix ---\n');
disp('Ob = [C; CA] ='); disp(Ob);
fprintf('Rank = %d -> %s\n', rank(Ob), iif(rank(Ob)==n, 'OBSERVABLE', 'NOT observable'));

%% Method 2 – Manual Observability Matrix
CA = C * A;
Ob_hand = [C; CA];
fprintf('\n--- Method 2: Manual Observability Matrix ---\n');
disp('Ob_hand = [C; CA] ='); disp(Ob_hand);
fprintf('Rank = %d -> %s\n', rank(Ob_hand), iif(rank(Ob_hand)==n, 'OBSERVABLE', 'NOT observable'));

%% Method 3 – Jordan Form Test (Observability)
fprintf('\n--- Method 3: Jordan Form Test (Observability) ---\n');
disp('J = '); disp(J);
disp('C_j = '); disp(C_jordan);
fprintf('All columns of C_j corresponding to each Jordan block must be non-zero.\n');

%% Method 4 – Observability Gramian
fprintf('\n--- Method 4: Observability Gramian ---\n');
if all(real(eig(A)) < 0)
    Wo = gram(sys_lin, 'o');
    disp('Wo = '); disp(Wo);
    disp('eig(Wo) = '); disp(eig(Wo));
    if all(eig(Wo) > 0)
        fprintf('Positive definite -> OBSERVABLE\n');
    else
        fprintf('Not positive definite -> NOT observable\n');
    end
else
    fprintf('System unstable, Gramian not defined.\n');
end

%% Method 5 – PBH Test for Observability
fprintf('\n--- Method 5: PBH Test (Observability) ---\n');
for i = 1:length(eigvals)
    lambda = eigvals(i);
    H_obs = [A - lambda*eye(n); C];
    fprintf('λ = %.4f + %.4fi, rank([A-λI; C]) = %d\n', ...
            real(lambda), imag(lambda), rank(H_obs));
end
fprintf('All ranks must be n (= %d) for observability.\n', n);

%% ===========================================================================
%  SECTION 10: Model Validation – Nonlinear vs Linear Response
% ===========================================================================
fprintf('\n========== 10. MODEL VALIDATION ==========\n');

% Simulation settings
dt = 0.001;
t_final = 20;
t = 0:dt:t_final;
x0 = x_ss;                    % start from equilibrium

% Nonlinear model as an anonymous function
f_nl = @(x, u) [...
    (x60_nom/x1)*(x20 - x(1)) - c0*x(1)*exp(-1/x(2));
    (x60_nom/x1)*(x30 - x(2)) + c1*x(1)*exp(-1/x(2)) - ...
        (c2*c3*u*(x(2) - x40)) / (x1*(c3*u + c4)) ];

% --- Test 1: Step Response ---
fprintf('--- Test 1: Step Input ---\n');
u_step = u_nom * ones(size(t));
step_mag = 0.1;
u_step(t >= 1) = u_nom + step_mag;   % 10% step at τ = 1

% Simulate nonlinear
x_nl_step = zeros(2, length(t));
x_nl_step(:,1) = x0;
for k = 1:length(t)-1
    dx = f_nl(x_nl_step(:,k), u_step(k));
    x_nl_step(:,k+1) = x_nl_step(:,k) + dx * dt;
end
y_nl_step = x_nl_step(2,:);

% Simulate linear
y_lin_step = lsim(sys_lin, u_step - u_nom, t, zeros(size(x0))) + x_ss(2);
y_lin_step = y_lin_step(:)';

rmse_step = sqrt(mean((y_nl_step - y_lin_step).^2));
fprintf('RMSE (step) = %.6f\n', rmse_step);

% --- Test 2: Sinusoidal Response ---
fprintf('--- Test 2: Sinusoidal Input ---\n');
freq = 0.5;          % rad/s
amp_sin = 0.05;
u_sin = u_nom + amp_sin * sin(freq * t);

% Simulate nonlinear
x_nl_sin = zeros(2, length(t));
x_nl_sin(:,1) = x0;
for k = 1:length(t)-1
    dx = f_nl(x_nl_sin(:,k), u_sin(k));
    x_nl_sin(:,k+1) = x_nl_sin(:,k) + dx * dt;
end
y_nl_sin = x_nl_sin(2,:);

% Simulate linear
y_lin_sin = lsim(sys_lin, u_sin - u_nom, t, zeros(size(x0))) + x_ss(2);
y_lin_sin = y_lin_sin(:)';

rmse_sin = sqrt(mean((y_nl_sin - y_lin_sin).^2));
fprintf('RMSE (sinusoidal) = %.6f\n', rmse_sin);

% --- Plot comparison ---
figure;
subplot(2,1,1);
plot(t, y_nl_step, 'b-', 'LineWidth', 1.5); hold on;
plot(t, y_lin_step, 'r--', 'LineWidth', 1.5);
legend('Nonlinear', 'Linear', 'Location', 'best');
xlabel('\tau'); ylabel('x_3');
title(sprintf('Step Response Comparison (RMSE = %.2e)', rmse_step));
grid on;

subplot(2,1,2);
plot(t, y_nl_sin, 'b-', 'LineWidth', 1.5); hold on;
plot(t, y_lin_sin, 'r--', 'LineWidth', 1.5);
legend('Nonlinear', 'Linear', 'Location', 'best');
xlabel('\tau'); ylabel('x_3');
title(sprintf('Sinusoidal Response (f = %.1f rad/s, RMSE = %.2e)', freq, rmse_sin));
grid on;

%% ===========================================================================
%  SECTION 11: RMSE Analysis for Varying Input Parameters
% ===========================================================================
fprintf('\n========== 11. RMSE VS INPUT AMPLITUDE & FREQUENCY ==========\n');

% Reuse dt, t, x0, f_nl, sys_lin (ensure they exist)
if ~exist('dt','var'), dt = 0.01; t_final = 20; t = 0:dt:t_final; end
if ~exist('x0','var'), x0 = x_ss; end
if ~exist('sys_lin','var'), sys_lin = ss(A, B, C, D); end

% --- 1. Step input: vary amplitude ---
step_amps = 0.05:0.05:1;
rmse_step_amps = zeros(size(step_amps));
for i = 1:length(step_amps)
    amp = step_amps(i);
    u_step = u_nom * ones(size(t));
    u_step(t >= 1) = u_nom + amp;
    % Nonlinear
    x_nl = zeros(2, length(t)); x_nl(:,1) = x0;
    for k = 1:length(t)-1
        dx = f_nl(x_nl(:,k), u_step(k));
        x_nl(:,k+1) = x_nl(:,k) + dx * dt;
    end
    y_nl = x_nl(2,:);
    % Linear
    y_lin = lsim(sys_lin, u_step - u_nom, t, zeros(size(x0))) + x_ss(2);
    y_lin = y_lin(:)';
    rmse_step_amps(i) = sqrt(mean((y_nl - y_lin).^2));
end

% --- 2. Sinusoidal input: vary amplitude (fixed freq = 1 Hz) ---
sin_freq_Hz = 1;
sin_freq_rad = 2*pi*sin_freq_Hz;
sin_amps = 0.05:0.05:1;
rmse_sin_amps = zeros(size(sin_amps));
for i = 1:length(sin_amps)
    amp = sin_amps(i);
    u_sin = u_nom + amp * sin(sin_freq_rad * t);
    x_nl = zeros(2, length(t)); x_nl(:,1) = x0;
    for k = 1:length(t)-1
        dx = f_nl(x_nl(:,k), u_sin(k));
        x_nl(:,k+1) = x_nl(:,k) + dx * dt;
    end
    y_nl = x_nl(2,:);
    y_lin = lsim(sys_lin, u_sin - u_nom, t, zeros(size(x0))) + x_ss(2);
    y_lin = y_lin(:)';
    rmse_sin_amps(i) = sqrt(mean((y_nl - y_lin).^2));
end

% --- 3. Sinusoidal input: vary frequency (fixed amp = 0.1) ---
fixed_amp = 0.1;
freqs_Hz = 0.5:0.5:50;
rmse_sin_freqs = zeros(size(freqs_Hz));
for i = 1:length(freqs_Hz)
    freq_rad = 2*pi*freqs_Hz(i);
    u_sin = u_nom + fixed_amp * sin(freq_rad * t);
    x_nl = zeros(2, length(t)); x_nl(:,1) = x0;
    for k = 1:length(t)-1
        dx = f_nl(x_nl(:,k), u_sin(k));
        x_nl(:,k+1) = x_nl(:,k) + dx * dt;
    end
    y_nl = x_nl(2,:);
    y_lin = lsim(sys_lin, u_sin - u_nom, t, zeros(size(x0))) + x_ss(2);
    y_lin = y_lin(:)';
    rmse_sin_freqs(i) = sqrt(mean((y_nl - y_lin).^2));
end

% --- Plot RMSE analyses ---
figure;

subplot(3,1,1);
plot(step_amps, rmse_step_amps, 'b-o', 'LineWidth', 1.5);
xlabel('Step Amplitude \Deltau'); ylabel('RMSE');
title('RMSE vs Step Input Amplitude'); grid on;

subplot(3,1,2);
plot(sin_amps, rmse_sin_amps, 'r-o', 'LineWidth', 1.5);
xlabel('Sinusoidal Amplitude \Deltau'); ylabel('RMSE');
title(sprintf('RMSE vs Sinusoidal Amplitude (f = %.1f Hz)', sin_freq_Hz)); grid on;

subplot(3,1,3);
semilogy(freqs_Hz, rmse_sin_freqs, 'g-o', 'LineWidth', 1.5);
xlabel('Frequency (Hz)'); ylabel('RMSE');
title(sprintf('RMSE vs Sinusoidal Frequency (A = %.1f)', fixed_amp)); grid on;

%% ===========================================================================
%  Helper inline function for conditional strings
% ===========================================================================
function s = iif(condition, str_true, str_false)
    if condition
        s = str_true;
    else
        s = str_false;
    end
end