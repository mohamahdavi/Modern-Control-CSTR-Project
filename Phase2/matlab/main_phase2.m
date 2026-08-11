%% ===========================================================================
%  Phase 2: Control System Design, Observer, Robustness & Saturation
%  Non-isothermal Continuous Stirred Tank Reactor (CSTR)
%  Based on the linearized model from Phase 1
% ===========================================================================
%  This script performs:
%    1. State feedback control design and domain of attraction
%    2. Canonical form realizations (controller, observer, Jordan, etc.)
%    3. Servo control with integral action (tracking & disturbance rejection)
%    4. Observer design (Luenberger & integral observer, speed analysis)
%    5. Observer-based control vs. ideal state feedback comparison
%    6. Input saturation effects with pole scaling and anti-windup
%    7. Robustness analysis under parameter uncertainty (±30% in F_os)
% ===========================================================================
clear; clc; close all;

%% ===========================================================================
%  Shared Parameters and Linearized Model (from Phase 1)
% ===========================================================================
fprintf('========== SHARED PARAMETERS & LINEARIZED MODEL ==========\n');

% Physical constants (same as Phase 1)
V_s    = 1.36;      % m^3
F_os   = 1.13;      % m^3/h
C_aor  = 3.92;      % kmol/m^3
C_ao   = 8;         % kmol/m^3
alpha  = 7.08e10;   % 1/h
H_rxn  = 69815;     % kJ/kmol (equal to E_a)
E_a    = 69815;     % kJ/kmol
R_gas  = 8.314;     % kJ/(kmol·K)
rho    = 800;       % kg/m^3
Cp     = 3.13;      % kJ/(kg·K)
U      = 3065;      % kJ/(h·m^2·K)
A_ht   = 23.22;     % m^2
F_js   = 1.4130;    % m^3/h
V_j    = 0.085;     % m^3
rho_j  = 1000;      % kg/m^3
Cp_j   = 4.18;      % kJ/(kg·K)
T_o    = 294.7;     % K
T_jo   = 294.7;     % K

% Dimensionless parameters (Table 2)
c0 = V_s * alpha / F_os;
c1 = (V_s * alpha * R_gas * C_aor) / (F_os * rho * Cp);
c2 = U * A_ht / (rho * Cp * F_os);
c3 = V_s * F_js / (F_os * V_j);
c4_star = (rho * Cp * V_s) / (rho_j * Cp_j * V_j);
c4 = c4_star * c2;

% Disturbances / constant terms
x1       = 1;            % V/V_s = 1
x20      = C_ao / C_aor; % dimensionless inlet concentration
x30      = R_gas * T_o / E_a;          % dimensionless inlet temperature
x40      = R_gas * T_jo / E_a;         % dimensionless inlet cooling water temp.
x60_nom  = 1;            % F_o / F_os = 1 (nominal)
u_nom    = 1;            % nominal dimensionless cooling flow

% Equilibrium point (from Phase 1)
x2_ss = 1.941355;
x3_ss = 0.035537;
x_ss = [x2_ss; x3_ss];

% Linearized state-space matrices (from Phase 1)
A = [-1.0512, -78.7580;
      0.0007,  -1.9024];
B = [0;
     -0.000787];
C = [0, 1];
D = 0;

fprintf('c0 = %.6f, c1 = %.6f, c2 = %.6f, c3 = %.6f, c4 = %.6f\n', c0,c1,c2,c3,c4);
fprintf('Equilibrium: x2 = %.6f, x3 = %.6f\n', x2_ss, x3_ss);
fprintf('Linearized system: A = \n'); disp(A);
fprintf('B = \n'); disp(B);
fprintf('C = \n'); disp(C); fprintf('\n');

%% ===========================================================================
%  SECTION 1: State Feedback Control and Domain of Attraction
% ===========================================================================
fprintf('========== 1. STATE FEEDBACK CONTROL ==========\n');

% Desired closed-loop poles
des_poles = [-3.0, -4.5];
K = place(A, B, des_poles);
fprintf('State feedback gain K = [%.4f, %.4f]\n', K(1), K(2));
fprintf('Closed-loop poles: %.2f, %.2f\n\n', des_poles(1), des_poles(2));

% Closed-loop nonlinear model for simulation
f_cl = @(t,x) cstr_sf(x, K, x_ss, u_nom, x1, x20, x30, x40, x60_nom, c0, c1, c2, c3, c4);

% Time responses from different initial conditions
tspan = [0 25];
X0 = {[1.95; 0.036], [2.30; 0.045], [2.80; 0.060], [1.20; 0.025], [3.20; 0.090]};
lbl = {'Near','Medium','Far','Opposite','Outside'};
cols = lines(5);

figure('Name','S1: Time Responses','Position',[100 100 900 600]);
h1a = gobjects(1,5); h1b = gobjects(1,5);
for i = 1:5
    [t,x] = ode45(f_cl, tspan, X0{i});
    subplot(2,1,1); hold on;
    h1a(i) = plot(t, x(:,1), 'Color',cols(i,:), 'LineWidth',1.5);
    subplot(2,1,2); hold on;
    h1b(i) = plot(t, x(:,2), 'Color',cols(i,:), 'LineWidth',1.5);
end
subplot(2,1,1); ylabel('x_2'); grid on; legend(h1a, lbl, 'Location','best'); title('Concentration');
subplot(2,1,2); xlabel('\tau'); ylabel('x_3'); grid on; legend(h1b, lbl, 'Location','best'); title('Temperature');

% Grid search for Domain of Attraction
fprintf('Estimating domain of attraction via grid search...\n');
x2g = 0.6:0.08:3.2; x3g = 0.015:0.003:0.10;
[X2,X3] = meshgrid(x2g, x3g);
Stable = zeros(size(X2));
for i = 1:numel(X2)
    try
        [~,x] = ode45(f_cl, [0 30], [X2(i); X3(i)]);
        if all(isfinite(x(end,:))) && norm(x(end,:)'-x_ss) < 0.05
            Stable(i) = 1;
        end
    end
end

figure('Name','S1: DoA','Position',[150 150 800 600]);
contourf(X2, X3, Stable, [1 1], 'LineColor','none');
colormap([1 0.7 0.7; 0.6 1 0.6]); hold on;
h_eq = plot(x2_ss, x3_ss, 'kp', 'MarkerSize',14, 'MarkerFaceColor','k');
xlabel('x_2'); ylabel('x_3'); title('Domain of Attraction');
legend(h_eq, {'Equilibrium'}, 'Location','best'); grid on; box on;

% Phase portrait
figure('Name','S1: Phase Portrait','Position',[200 200 800 600]); hold on; grid on;
h_traj = gobjects(1,5);
for i = 1:5
    [~,x] = ode45(f_cl, [0 25], X0{i});
    h_traj(i) = plot(x(:,1), x(:,2), 'LineWidth',1.8, 'Color',cols(i,:));
    plot(X0{i}(1), X0{i}(2), 'o', 'Color',cols(i,:), 'MarkerSize',8, ...
        'MarkerFaceColor',cols(i,:), 'HandleVisibility','off');
end
h_eq2 = plot(x2_ss, x3_ss, 'kp', 'MarkerSize',14, 'MarkerFaceColor','k');
xlabel('x_2'); ylabel('x_3'); title('Phase Portrait'); axis([0.5 3.3 0.01 0.11]);
legend([h_traj, h_eq2], [lbl, {'Equilibrium'}], 'Location','best');

%% ===========================================================================
%  SECTION 2: Canonical Form Realizations
% ===========================================================================
fprintf('\n========== 2. CANONICAL FORMS ==========\n');

% Transfer function
[num, den] = ss2tf(A, B, C, D);
n1 = num(2); n2 = num(3);
a1 = den(2); a2 = den(3);
fprintf('Transfer function: (%.4g s + %.4g) / (s^2 + %.4g s + %.4g)\n', n1, n2, a1, a2);

% --- Controller canonical form (bottom-companion) ---
A_ccf = [0 1; -a2 -a1];
B_ccf = [0; 1];
C_ccf = [n2 n1];
fprintf('\nController Canonical Form:\nA = '); disp(A_ccf);
fprintf('B = '); disp(B_ccf);
fprintf('C = '); disp(C_ccf);

% --- Observer canonical form (exact dual) ---
A_ocf = A_ccf.';
B_ocf = C_ccf.';
C_ocf = B_ccf.';
fprintf('\nObserver Canonical Form:\nA = '); disp(A_ocf);
fprintf('B = '); disp(B_ocf);
fprintf('C = '); disp(C_ocf);

% --- Alternative controller canonical form (top-companion) ---
A_ccf2 = [-a1 -a2; 1 0];
B_ccf2 = [1; 0];
C_ccf2 = [n1 n2];
fprintf('\nController Canonical Form (alt):\nA = '); disp(A_ccf2);
fprintf('B = '); disp(B_ccf2);
fprintf('C = '); disp(C_ccf2);

% --- Dual observer form of the alternative ---
A_ocf2 = A_ccf2.';
B_ocf2 = C_ccf2.';
C_ocf2 = B_ccf2.';
fprintf('\nObserver Canonical Form (alt):\nA = '); disp(A_ocf2);
fprintf('B = '); disp(B_ocf2);
fprintf('C = '); disp(C_ocf2);

% --- Jordan form via eigen-decomposition ---
[V, Dm] = eig(A);
A_jor = Dm;
B_jor = V \ B;
C_jor = C * V;
fprintf('\nJordan Form:\nJ = '); disp(A_jor);
fprintf('B_j = '); disp(B_jor);
fprintf('C_j = '); disp(C_jor);

% Collect all realizations
names = {'Original','Controller CF','Observer CF','Controller CF (alt)','Observer CF (alt)','Jordan'};
Ac = {A,      A_ccf,  A_ocf,  A_ccf2,  A_ocf2,  A_jor};
Bc = {B,      B_ccf,  B_ocf,  B_ccf2,  B_ocf2,  B_jor};
Cc = {C,      C_ccf,  C_ocf,  C_ccf2,  C_ocf2,  C_jor};
Dc = {D, 0, 0, 0, 0, 0};

% Verify transfer functions
fprintf('\nTransfer Function Verification:\n');
for i = 1:length(names)
    [numi,deni] = ss2tf(Ac{i},Bc{i},Cc{i},Dc{i});
    fprintf('%-20s num=[%+.3e %+.3e %+.3e] den=[1 %+.4f %+.4f]\n', ...
        names{i}, numi(1), numi(2), numi(3), deni(2), deni(3));
end

% Step responses (output and states)
n = length(names); cols2 = lines(n); lst = {'-','--','-.',':','-','--'};
t = 0:0.01:10; u = ones(size(t));
figure('Name','S2: Step Responses','Position',[100 100 1000 700]);
h2 = gobjects(1,n);
for i = 1:n
    [y,~,x] = lsim(ss(Ac{i},Bc{i},Cc{i},Dc{i}), u, t);
    subplot(3,1,1); hold on;
    h2(i) = plot(t,y,lst{i},'Color',cols2(i,:),'LineWidth',1.2);
    subplot(3,1,2); hold on; plot(t,x(:,1),lst{i},'Color',cols2(i,:),'LineWidth',1.2);
    subplot(3,1,3); hold on; plot(t,x(:,2),lst{i},'Color',cols2(i,:),'LineWidth',1.2);
end
subplot(3,1,1); ylabel('y'); title('Output (Identical)'); grid on;
legend(h2, names, 'Location','best');
subplot(3,1,2); ylabel('x_1'); title('State x_1 (Different)'); grid on;
subplot(3,1,3); xlabel('t'); ylabel('x_2'); title('State x_2 (Different)'); grid on;

%% ===========================================================================
%  SECTION 3: Servo Control with Integral Action
% ===========================================================================
fprintf('\n========== 3. SERVO CONTROL WITH INTEGRAL ACTION ==========\n');

% Augmented system: [A 0; -C 0]
A_aug = [A, zeros(2,1); -C, 0];
B_aug = [B; 0];
K_aug = place(A_aug, B_aug, [-4.0, -5.0, -1.2]);
K_s = K_aug(1:2); Ki = K_aug(3);
fprintf('Augmented state feedback: K_s = [%.4f, %.4f], Ki = %.4f\n', K_s(1), K_s(2), Ki);

% Reference and disturbance timing
delta_r = 0.0005;   % step change in reference
t_ref = 2;           % time of reference step
t_dist = 15;         % time of disturbance
d_amp = 0.10;        % disturbance amplitude
r_fun = @(tt) x3_ss + delta_r*(tt >= t_ref);
d_fun = @(tt) d_amp*(tt >= t_dist);

% Classical state feedback simulation
[t_c,x_c] = ode45(@(t,x) cstr_classical(t,x,K,x_ss,u_nom,d_fun,x1,x20,x30,x40,x60_nom,c0,c1,c2,c3,c4), [0 40], x_ss);
y_c = x_c(:,2);
r_c = arrayfun(r_fun,t_c);
e_c = r_c' - y_c;

% Servo (integral) simulation
[t_s,x_s] = ode45(@(t,x) cstr_servo(t,x,K_s,Ki,x_ss,u_nom,r_fun,d_fun,x1,x20,x30,x40,x60_nom,c0,c1,c2,c3,c4), [0 40], [x_ss; 0]);
y_s = x_s(:,2);
r_s = arrayfun(r_fun,t_s);
e_s = r_s' - y_s;

fprintf('Steady-state |error| (t>35s): Classical=%.2e  Servo=%.2e\n', ...
    mean(abs(e_c(t_c>35))), mean(abs(e_s(t_s>35))));

% Plot comparison
figure('Name','S3: Classical vs Servo','Position',[150 150 900 500]);
subplot(2,1,1); hold on; grid on;
plot(t_c, y_c, 'r-', 'LineWidth',1.5);
plot(t_s, y_s, 'b-', 'LineWidth',1.5);
plot(t_c, r_c, 'k--', 'LineWidth',1.2);
xline(t_dist,'m:','Disturbance','HandleVisibility','off');
ylabel('x_3'); title('Output');
legend({'Classical','Servo','Reference'}, 'Location','best');

subplot(2,1,2); hold on; grid on;
plot(t_c, e_c, 'r-', 'LineWidth',1.5);
plot(t_s, e_s, 'b-', 'LineWidth',1.5);
yline(0, 'k--','HandleVisibility','off');
xline(t_dist,'m:','Disturbance','HandleVisibility','off');
xlabel('t'); ylabel('e(t)'); title('Tracking Error');

figure('Name','S3: Servo Signals','Position',[200 200 900 400]); hold on; grid on;
plot(t_s, x_s(:,3), 'c--', 'LineWidth',1.2);
xline(t_ref,'g:','Reference','HandleVisibility','off');
xline(t_dist,'r:','Disturbance','HandleVisibility','off');
xlabel('t'); ylabel('\xi(t)'); title('Integrator State');
legend({'Integrator \xi'}, 'Location','best');

%% ===========================================================================
%  SECTION 4: Observer Design and State Estimation
% ===========================================================================
fprintf('\n========== 4. OBSERVER DESIGN ==========\n');

% Simulation settings for all observer tests
dt = 1e-4;
T_end = 20;
t = 0:dt:T_end;
steps = length(t);

% --- 4.1 Luenberger observer (no integral) ---
observer_poles = [-9, -12];
L = place(A', C', observer_poles)';
fprintf('\n--- 4.1 Luenberger Observer (no integral) ---\n');
fprintf('Desired observer poles: [%.0f, %.0f]\n', observer_poles(1), observer_poles(2));
fprintf('Observer gain L = \n'); disp(L);

% Simulation without integral
x_real = zeros(steps, 2); x_hat = zeros(steps, 2);
x_real(1,:) = x_ss'; x_hat(1,:) = x_ss';

for k = 1:steps-1
    x2 = x_real(k,1); x3 = x_real(k,2);
    x_hat_k = x_hat(k,:)';
    
    d_val = 0.05 * (t(k) >= 5);
    
    % Control using estimated state
    u_hat = u_nom - K(1)*(x_hat_k(1)-x_ss(1)) - K(2)*(x_hat_k(2)-x_ss(2));
    u_real = u_hat + d_val;
    
    % Real nonlinear dynamics
    dx2_real = (x60_nom/x1)*(x20-x2) - c0*x2*exp(-1/x3);
    dx3_real = (x60_nom/x1)*(x30-x3) + c1*x2*exp(-1/x3) - (c2*c3*u_real*(x3-x40))/(x1*(c3*u_real+c4));
    
    y_meas = C * [x2; x3];
    error_y = y_meas - C * x_hat_k;
    
    % Observer (nonlinear model + correction)
    xh1 = x_hat_k(1); xh2 = x_hat_k(2);
    dxh1 = (x60_nom/x1)*(x20 - xh1) - c0*xh1*exp(-1/xh2);
    dxh2 = (x60_nom/x1)*(x30 - xh2) + c1*xh1*exp(-1/xh2) - (c2*c3*u_hat*(xh2 - x40))/(x1*(c3*u_hat + c4));
    
    dx_hat = [dxh1; dxh2] + L * error_y;
    
    x_real(k+1,:) = x_real(k,:) + dt * [dx2_real, dx3_real];
    x_hat(k+1,:) = x_hat(k,:) + dt * dx_hat';
end
fprintf('Estimation error (final): x1_err=%.4e, x2_err=%.4e\n', ...
    x_real(end,1)-x_hat(end,1), x_real(end,2)-x_hat(end,2));

% --- 4.2 Integral Observer (disturbance rejection) ---
% Place poles for augmented observer system
Lp = [0; 19.0464];
Li = [-6282.869; 106.923];
A_obs_int = [A - Lp*C, Li; -C, 0];
fprintf('\n--- 4.2 Integral Observer ---\n');
fprintf('Lp = \n'); disp(Lp);
fprintf('Li = \n'); disp(Li);
fprintf('Observer eigenvalues: \n'); disp(eig(A_obs_int));

% Simulation with integral observer
x_real_int = zeros(steps, 2); x_hat_int = zeros(steps, 2); z = zeros(steps,1);
x_real_int(1,:) = x_ss'; x_hat_int(1,:) = x_ss';

for k = 1:steps-1
    x2 = x_real_int(k,1); x3 = x_real_int(k,2);
    x_hat_k = x_hat_int(k,:)';
    
    d_val = 0.05 * (t(k) >= 5);
    u_hat = u_nom - K(1)*(x_hat_k(1)-x_ss(1)) - K(2)*(x_hat_k(2)-x_ss(2));
    u_real = u_hat + d_val;
    
    dx2_real = (x60_nom/x1)*(x20-x2) - c0*x2*exp(-1/x3);
    dx3_real = (x60_nom/x1)*(x30-x3) + c1*x2*exp(-1/x3) - (c2*c3*u_real*(x3-x40))/(x1*(c3*u_real+c4));
    
    y_meas = C * [x2; x3];
    error_y = y_meas - C * x_hat_k;
    
    xh1 = x_hat_k(1); xh2 = x_hat_k(2);
    dxh1 = (x60_nom/x1)*(x20 - xh1) - c0*xh1*exp(-1/xh2);
    dxh2 = (x60_nom/x1)*(x30 - xh2) + c1*xh1*exp(-1/xh2) - (c2*c3*u_real*(xh2 - x40))/(x1*(c3*u_real + c4));
    
    dx_hat = [dxh1; dxh2] + Lp * error_y + Li * z(k);
    dz = error_y;
    
    x_real_int(k+1,:) = x_real_int(k,:) + dt * [dx2_real, dx3_real];
    x_hat_int(k+1,:) = x_hat_int(k,:) + dt * dx_hat';
    z(k+1) = z(k) + dt * dz;
end
fprintf('Integral observer estimation error (final): x1_err=%.4e, x2_err=%.4e\n', ...
    x_real_int(end,1)-x_hat_int(end,1), x_real_int(end,2)-x_hat_int(end,2));

% Plot comparison of no-integral vs integral observer
figure;
subplot(2,2,1);
plot(t, x_real(:,1), 'b-', 'LineWidth',1.5); hold on;
plot(t, x_hat(:,1), 'r-', 'LineWidth',1.5);
xline(5, 'k:'); ylabel('x_1'); grid on; legend('Real','Estimated'); title('State 1 (no integral)');
subplot(2,2,2);
plot(t, x_real(:,2), 'b-', 'LineWidth',1.5); hold on;
plot(t, x_hat(:,2), 'r-', 'LineWidth',1.5);
xline(5, 'k:'); ylabel('x_2'); grid on; legend('Real','Estimated'); title('State 2 (no integral)');
subplot(2,2,3);
plot(t, x_real_int(:,1), 'b-', 'LineWidth',1.5); hold on;
plot(t, x_hat_int(:,1), 'r-', 'LineWidth',1.5);
xline(5, 'k:'); ylabel('x_1'); grid on; legend('Real','Estimated'); title('State 1 (with integral)');
subplot(2,2,4);
plot(t, x_real_int(:,2), 'b-', 'LineWidth',1.5); hold on;
plot(t, x_hat_int(:,2), 'r-', 'LineWidth',1.5);
xline(5, 'k:'); ylabel('x_2'); grid on; legend('Real','Estimated'); title('State 2 (with integral)');

% --- 4.3 Observer Speed Analysis (with noise) ---
noise_power = 1e-12;

% Gains for different observer speeds (integral observer)
Lp_slow  = [0; -0.7336];        Li_slow  = [113.31; 0.02728];
Lp_mid   = [5976.8918; 19.0464]; Li_mid   = [9050; 100];
Lp_fast  = [0; 217.0464];        Li_fast  = [-3.576e6; 12669.785];
Lp_slowest = [0; -2.7336];       Li_slowest = [-1248.61; 0.8356];

fprintf('\n--- 4.3 Observer Speed Comparison (with measurement noise) ---\n');
fprintf('Slow observer: Lp = [%.4f; %.4f], Li = [%.4f; %.4f]\n', Lp_slow(1),Lp_slow(2),Li_slow(1),Li_slow(2));
fprintf('Mid observer:  Lp = [%.4f; %.4f], Li = [%.4f; %.4f]\n', Lp_mid(1),Lp_mid(2),Li_mid(1),Li_mid(2));
fprintf('Fast observer: Lp = [%.4f; %.4f], Li = [%.4f; %.4f]\n', Lp_fast(1),Lp_fast(2),Li_fast(1),Li_fast(2));

% Anonymous function wrapper (returns 5 outputs)
simulate_observer = @(Lp,Li,label) simulate_integral_observer(...
    steps, dt, t, x_ss, A, B, C, u_nom, K, noise_power, ...
    x1,x20,x30,x40,x60_nom,c0,c1,c2,c3,c4, Lp, Li);

% Run all observer speed simulations
[x_slow, xh_slow, ~, u_slow, y_slow] = simulate_observer(Lp_slow, Li_slow, 'slow');
[x_mid,  xh_mid,  ~, u_mid,  y_mid]  = simulate_observer(Lp_mid,  Li_mid,  'mid');
[x_fast, xh_fast, ~, u_fast, y_fast] = simulate_observer(Lp_fast, Li_fast, 'fast');
[x_slowest, xh_slowest, ~, u_slowest, y_slowest] = simulate_observer(Lp_slowest, Li_slowest, 'slowest');

% Plot all speed results in 3x4 grid (slow, mid, fast)
figure;
plot_observer_speed(t, u_slow, y_slow, x_slow, xh_slow, '10x slower', 1);
plot_observer_speed(t, u_mid,  y_mid,  x_mid,  xh_mid,  'mid', 5);
plot_observer_speed(t, u_fast, y_fast, x_fast, xh_fast, '10x faster', 9);

% For "slowest" and "slow" compare in a separate 2x4 grid
figure;
subplot(2,4,1); plot(t,u_slowest,'g-'); xlabel('Time'); ylabel('u'); title('Control (100x slower)'); grid on; xlim([0 19.99]);
subplot(2,4,2); plot(t,y_slowest,'r-'); xlabel('Time'); ylabel('y'); title('Output (100x slower)'); grid on; xlim([0 19.99]);
subplot(2,4,3); plot(t,x_slowest(:,1),'b-',t,xh_slowest(:,1),'r--'); xlabel('Time'); ylabel('x_1'); title('State 1 (100x slower)'); legend('Real','Est'); grid on; xlim([0 19.99]);
subplot(2,4,4); plot(t,x_slowest(:,2),'b-',t,xh_slowest(:,2),'r--'); xlabel('Time'); ylabel('x_2'); title('State 2 (100x slower)'); legend('Real','Est'); grid on; xlim([0 19.99]);
subplot(2,4,5); plot(t,u_slow,'g-'); xlabel('Time'); ylabel('u'); title('Control (10x slower)'); grid on; xlim([0 19.99]);
subplot(2,4,6); plot(t,y_slow,'r-'); xlabel('Time'); ylabel('y'); title('Output (10x slower)'); grid on; xlim([0 19.99]);
subplot(2,4,7); plot(t,x_slow(:,1),'b-',t,xh_slow(:,1),'r--'); xlabel('Time'); ylabel('x_1'); title('State 1 (10x slower)'); legend('Real','Est'); grid on; xlim([0 19.99]);
subplot(2,4,8); plot(t,x_slow(:,2),'b-',t,xh_slow(:,2),'r--'); xlabel('Time'); ylabel('x_2'); title('State 2 (10x slower)'); legend('Real','Est'); grid on; xlim([0 19.99]);

%% ===========================================================================
%  SECTION 5: Observer-Based Control vs. Ideal State Feedback
% ===========================================================================
fprintf('\n========== 5. OBSERVER-BASED CONTROL VS. IDEAL ==========\n');

% Ideal state feedback (full state known)
x_ideal = zeros(steps,2); u_ideal = zeros(steps,1); y_ideal = zeros(steps,1);
x_ideal(1,:) = x_ss';
for k = 1:steps-1
    x2 = x_ideal(k,1); x3 = x_ideal(k,2);
    u = u_nom - K(1)*(x2-x_ss(1)) - K(2)*(x3-x_ss(2));
    u_ideal(k) = u;
    dx2 = (x60_nom/x1)*(x20-x2) - c0*x2*exp(-1/x3);
    dx3 = (x60_nom/x1)*(x30-x3) + c1*x2*exp(-1/x3) - (c2*c3*u*(x3-x40))/(x1*(c3*u+c4));
    x_ideal(k+1,:) = x_ideal(k,:) + dt * [dx2, dx3];
    y_ideal(k) = x3;
end
y_ideal(end) = x_ideal(end,2); u_ideal(end) = u_ideal(end-1);

% Observer-based control (using slow observer from Section 4)
x_obs = zeros(steps,2); x_hat_obs = zeros(steps,2); z_obs = zeros(steps,1);
u_obs = zeros(steps,1); y_obs = zeros(steps,1);
x_obs(1,:) = x_ss'; x_hat_obs(1,:) = x_ss' + [0.001, 0.0004]; % initial error
for k = 1:steps-1
    x2 = x_obs(k,1); x3 = x_obs(k,2);
    x_hat_k = x_hat_obs(k,:)';
    
    u_hat = u_nom - K(1)*(x_hat_k(1)-x_ss(1)) - K(2)*(x_hat_k(2)-x_ss(2));
    u_obs(k) = u_hat;
    
    dx2_real = (x60_nom/x1)*(x20-x2) - c0*x2*exp(-1/x3);
    dx3_real = (x60_nom/x1)*(x30-x3) + c1*x2*exp(-1/x3) - (c2*c3*u_hat*(x3-x40))/(x1*(c3*u_hat+c4));
    
    y_meas = C * [x2; x3]; y_obs(k) = y_meas;
    error_y = y_meas - C * x_hat_k;
    
    xh1 = x_hat_k(1); xh2 = x_hat_k(2);
    dxh1 = (x60_nom/x1)*(x20-xh1) - c0*xh1*exp(-1/xh2);
    dxh2 = (x60_nom/x1)*(x30-xh2) + c1*xh1*exp(-1/xh2) - (c2*c3*u_hat*(xh2-x40))/(x1*(c3*u_hat+c4));
    
    dx_hat = [dxh1; dxh2] + Lp_slow * error_y + Li_slow * z_obs(k);
    dz = error_y;
    
    x_obs(k+1,:) = x_obs(k,:) + dt * [dx2_real, dx3_real];
    x_hat_obs(k+1,:) = x_hat_obs(k,:) + dt * dx_hat';
    z_obs(k+1) = z_obs(k) + dt * dz;
end
y_obs(end) = x_obs(end,2); u_obs(end) = u_obs(end-1);

% Compare control and output
figure;
subplot(2,2,1);
plot(t, u_ideal, 'b-', 'LineWidth',1.5); xlabel('Time'); ylabel('u(t)'); title('Control Signal (Ideal)'); grid on; xlim([0 19.99]);
subplot(2,2,2);
plot(t, y_ideal, 'b-', 'LineWidth',1.5); xlabel('Time'); ylabel('y(t)'); title('Output (Ideal)'); grid on; xlim([0 19.99]);
subplot(2,2,3);
plot(t, u_obs, 'r-', 'LineWidth',1.5); xlabel('Time'); ylabel('u(t)'); title('Control Signal (Observer-based)'); grid on; xlim([0 19.99]);
subplot(2,2,4);
plot(t, y_obs, 'r-', 'LineWidth',1.5); xlabel('Time'); ylabel('y(t)'); title('Output (Observer-based)'); grid on; xlim([0 19.99]);

%% ===========================================================================
%  SECTION 6: Input Saturation and Anti-Windup
% ===========================================================================
fprintf('\n========== 6. INPUT SATURATION & ANTI-WINDUP ==========\n');

u_min = 0; u_max = 2;
scale_factors = [1, 1.5, 2, 2.5, 3, 3.5, 4, 5];
K_base = K; Lp_base = Lp_slow; Li_base = Li_slow;

% --- 6.1 Saturation only (no anti-windup) ---
figure('Name','Part 6: Control with Saturation');
for i = 1:length(scale_factors)
    sf = scale_factors(i);
    K_test = place(A, B, sf * [-3.0, -4.5]);
    Lp_test = sf * Lp_base;
    Li_test = sf * Li_base;
    
    % Simulate with saturation (returns 8 outputs)
    [~, ~, ~, u_sat, u_unsat, y_sat, ~, ~] = simulate_saturation(...
        steps, dt, t, x_ss, A, B, C, u_nom, u_min, u_max, ...
        K_test, Lp_test, Li_test, x1,x20,x30,x40,x60_nom,c0,c1,c2,c3,c4, false);
    
    subplot(4, 4, 2*(i-1)+1);
    plot(t, u_sat, 'g-', 'LineWidth',1.2); hold on;
    plot(t, u_unsat, 'b--', 'LineWidth',1);
    xlabel('Time'); ylabel('u(t)'); grid on;
    title(sprintf('Scale %g: u (saturation only)', sf));
    legend('Sat','Unsat','Location','best'); xlim([0 20]);
    
    subplot(4, 4, 2*(i-1)+2);
    plot(t, y_sat, 'r-', 'LineWidth',1.2);
    xlabel('Time'); ylabel('y(t)'); grid on;
    title(sprintf('Scale %g: y', sf));
    xlim([0 20]);
end

% --- 6.2 Anti-windup (freeze integrator when saturated) ---
figure('Name','Part 6: Control with Anti-Windup');
for i = 1:length(scale_factors)
    sf = scale_factors(i);
    K_test = place(A, B, sf * [-3.0, -4.5]);
    Lp_test = sf * Lp_base;
    Li_test = sf * Li_base;
    
    [~, ~, ~, u_aw, u_unsat_aw, y_aw, ~, ~] = simulate_saturation(...
        steps, dt, t, x_ss, A, B, C, u_nom, u_min, u_max, ...
        K_test, Lp_test, Li_test, x1,x20,x30,x40,x60_nom,c0,c1,c2,c3,c4, true);
    
    subplot(4, 4, 2*(i-1)+1);
    plot(t, u_aw, 'g-', 'LineWidth',1.2); hold on;
    plot(t, u_unsat_aw, 'b--', 'LineWidth',1);
    xlabel('Time'); ylabel('u(t)'); grid on;
    title(sprintf('Scale %g: u (AW)', sf));
    legend('Sat(AW)','Unsat','Location','best'); xlim([0 20]);
    
    subplot(4, 4, 2*(i-1)+2);
    plot(t, y_aw, 'r-', 'LineWidth',1.2);
    xlabel('Time'); ylabel('y(t)'); grid on;
    title(sprintf('Scale %g: y (AW)', sf));
    xlim([0 20]);
end

%% ===========================================================================
%  SECTION 7: Robustness to Parameter Uncertainty (F_os +30%)
% ===========================================================================
fprintf('\n========== 7. ROBUSTNESS ANALYSIS (F_os +30%%) ==========\n');

F_os_new = 1.3 * F_os;
c0_new = V_s*alpha/F_os_new;
c1_new = (V_s*alpha*R_gas*C_aor)/(F_os_new*rho*Cp);
c2_new = U*A_ht/(rho*Cp*F_os_new);
c3_new = V_s*F_js/(F_os_new*V_j);
c4_new = (rho*Cp*V_s)/(rho_j*Cp_j*V_j)*c2_new;

sf = 2;
K_test = place(A, B, sf * [-3.0, -4.5]);
Lp_test = sf * Lp_base;
Li_test = sf * Li_base;

% --- 7.1 Without saturation ---
fprintf('Simulating 30%% increase in F_os (no saturation)...\n');
[x_rob, xh_rob, u_rob] = simulate_robustness(steps, dt, t, x_ss, A, B, C, u_nom, ...
    K_test, Lp_test, Li_test, x1,x20,x30,x40,x60_nom,c0,c1,c2,c3,c4, ...
    c0_new,c1_new,c2_new,c3_new,c4_new, false);
fprintf('Final output estimation error: %.4e\n', x_rob(end,2)-xh_rob(end,2));

figure;
subplot(2,1,1); plot(t, u_rob, 'g-'); xlabel('Time'); ylabel('u'); title('Control (no saturation)'); grid on; xlim([0 20]);
subplot(2,1,2); plot(t, x_rob(:,2), 'b-', t, xh_rob(:,2), 'r--'); xlabel('Time'); ylabel('x_3'); title('Output'); grid on; xlim([0 20]); legend('Real','Est');

% --- 7.2 With saturation (same gains) ---
fprintf('Simulating with saturation (u in [0,2])...\n');
[x_rob2, xh_rob2, u_rob2] = simulate_robustness(steps, dt, t, x_ss, A, B, C, u_nom, ...
    K_test, Lp_test, Li_test, x1,x20,x30,x40,x60_nom,c0,c1,c2,c3,c4, ...
    c0_new,c1_new,c2_new,c3_new,c4_new, true);
fprintf('Final output estimation error (with saturation): %.4e\n', x_rob2(end,2)-xh_rob2(end,2));

figure;
subplot(2,1,1); plot(t, u_rob2, 'g-'); xlabel('Time'); ylabel('u'); title('Control (with saturation)'); grid on; xlim([0 20]);
subplot(2,1,2); plot(t, x_rob2(:,2), 'b-', t, xh_rob2(:,2), 'r--'); xlabel('Time'); ylabel('x_3'); title('Output'); grid on; xlim([0 20]); legend('Real','Est');

fprintf('\nPhase 2 analysis complete.\n');

%% ===========================================================================
%  Helper Functions
% ===========================================================================

function dx = cstr_sf(x, K, x_ss, u_nom, x1, x20, x30, x40, x60, c0, c1, c2, c3, c4)
    x2 = x(1); x3 = x(2);
    u = max(u_nom - K(1)*(x2-x_ss(1)) - K(2)*(x3-x_ss(2)), 0);
    dx2 = (x60/x1)*(x20-x2) - c0*x2*exp(-1/x3);
    dx3 = (x60/x1)*(x30-x3) + c1*x2*exp(-1/x3) - (c2*c3*u*(x3-x40))/(x1*(c3*u+c4));
    dx = [dx2; dx3];
end

function dx = cstr_classical(t, x, K, x_ss, u_nom, d_fun, x1, x20, x30, x40, x60, c0, c1, c2, c3, c4)
    x2 = x(1); x3 = x(2);
    u = max(u_nom - K(1)*(x2-x_ss(1)) - K(2)*(x3-x_ss(2)), 0);
    u_p = u + d_fun(t);
    dx2 = (x60/x1)*(x20-x2) - c0*x2*exp(-1/x3);
    dx3 = (x60/x1)*(x30-x3) + c1*x2*exp(-1/x3) - (c2*c3*u_p*(x3-x40))/(x1*(c3*u_p+c4));
    dx = [dx2; dx3];
end

function dx = cstr_servo(t, x_vec, K, Ki, x_ss, u_nom, r_fun, d_fun, x1, x20, x30, x40, x60, c0, c1, c2, c3, c4)
    x2 = x_vec(1); x3 = x_vec(2); xi = x_vec(3);
    u = max(u_nom - K(1)*(x2-x_ss(1)) - K(2)*(x3-x_ss(2)) - Ki*xi, 0);
    u_p = u + d_fun(t);
    dx2 = (x60/x1)*(x20-x2) - c0*x2*exp(-1/x3);
    dx3 = (x60/x1)*(x30-x3) + c1*x2*exp(-1/x3) - (c2*c3*u_p*(x3-x40))/(x1*(c3*u_p+c4));
    dxi = r_fun(t) - x3;
    dx = [dx2; dx3; dxi];
end

function [x_real, x_hat, z, u, u_unsat, y, xh_out, z_out] = simulate_saturation(steps, dt, t, x_ss, A, B, C, u_nom, u_min, u_max, K, Lp, Li, x1,x20,x30,x40,x60,c0,c1,c2,c3,c4, antiwindup)
    x_real = zeros(steps,2); x_hat = zeros(steps,2); z = zeros(steps,1);
    u = zeros(steps,1); u_unsat = zeros(steps,1); y = zeros(steps,1);
    x_real(1,:) = x_ss'; x_hat(1,:) = x_ss';
    for k = 1:steps-1
        x2 = x_real(k,1); x3 = x_real(k,2);
        x_hat_k = x_hat(k,:)';
        
        d_val = 0.05 * (t(k) >= 5);
        u_hat = u_nom - K(1)*(x_hat_k(1)-x_ss(1)) - K(2)*(x_hat_k(2)-x_ss(2));
        u_unsat(k) = u_hat;
        if u_hat > u_max, u_hat = u_max; elseif u_hat < u_min, u_hat = u_min; end
        u(k) = u_hat;
        u_real = u_hat + d_val;
        
        dx2_real = (x60/x1)*(x20-x2) - c0*x2*exp(-1/x3);
        dx3_real = (x60/x1)*(x30-x3) + c1*x2*exp(-1/x3) - (c2*c3*u_real*(x3-x40))/(x1*(c3*u_real+c4));
        
        y_meas = C * [x2; x3]; y(k) = y_meas;
        error_y = y_meas - C * x_hat_k;
        
        xh1 = x_hat_k(1); xh2 = x_hat_k(2);
        dxh1 = (x60/x1)*(x20-xh1) - c0*xh1*exp(-1/xh2);
        dxh2 = (x60/x1)*(x30-xh2) + c1*xh1*exp(-1/xh2) - (c2*c3*u_hat*(xh2-x40))/(x1*(c3*u_hat+c4));
        
        dx_hat = [dxh1; dxh2] + Lp * error_y + Li * z(k);
        dz = error_y;
        if antiwindup && (u_hat == u_max || u_hat == u_min), dz = 0; end
        
        x_real(k+1,:) = x_real(k,:) + dt * [dx2_real, dx3_real];
        x_hat(k+1,:) = x_hat(k,:) + dt * dx_hat';
        z(k+1) = z(k) + dt * dz;
    end
    y(end) = x_real(end,2); u(end) = u(end-1); u_unsat(end) = u_unsat(end-1);
    xh_out = x_hat; z_out = z;
end

function [x_real, x_hat, z, u, y] = simulate_integral_observer(steps, dt, t, x_ss, A, B, C, u_nom, K, noise_power, x1,x20,x30,x40,x60,c0,c1,c2,c3,c4, Lp, Li)
    x_real = zeros(steps,2); x_hat = zeros(steps,2); z = zeros(steps,1);
    u = zeros(steps,1); y = zeros(steps,1);
    x_real(1,:) = x_ss'; x_hat(1,:) = x_ss';
    for k = 1:steps-1
        x2 = x_real(k,1); x3 = x_real(k,2);
        x_hat_k = x_hat(k,:)';
        d_val = 0.05 * (t(k) >= 5);
        u_hat = u_nom - K(1)*(x_hat_k(1)-x_ss(1)) - K(2)*(x_hat_k(2)-x_ss(2));
        u_real = u_hat + d_val;
        u(k) = u_real;
        dx2_real = (x60/x1)*(x20-x2) - c0*x2*exp(-1/x3);
        dx3_real = (x60/x1)*(x30-x3) + c1*x2*exp(-1/x3) - (c2*c3*u_real*(x3-x40))/(x1*(c3*u_real+c4));
        y_meas = C * [x2; x3] + sqrt(noise_power)*randn;
        y(k) = y_meas;
        error_y = y_meas - C * x_hat_k;
        xh1 = x_hat_k(1); xh2 = x_hat_k(2);
        dxh1 = (x60/x1)*(x20-xh1) - c0*xh1*exp(-1/xh2);
        dxh2 = (x60/x1)*(x30-xh2) + c1*xh1*exp(-1/xh2) - (c2*c3*u_real*(xh2-x40))/(x1*(c3*u_real+c4));
        dx_hat = [dxh1; dxh2] + Lp*error_y + Li*z(k);
        x_real(k+1,:) = x_real(k,:) + dt*[dx2_real,dx3_real];
        x_hat(k+1,:) = x_hat(k,:) + dt*dx_hat';
        z(k+1) = z(k) + dt*error_y;
    end
    y(end) = x_real(end,2); u(end) = u(end-1);
end

function plot_observer_speed(t, u, y, x_real, x_hat, label, row_start)
    subplot(3,4,row_start);
    plot(t, u, 'g-'); xlabel('Time'); ylabel('u'); title(['Control (' label ')']); xlim([0 19.99]); grid on;
    subplot(3,4,row_start+1);
    plot(t, y, 'r-'); xlabel('Time'); ylabel('y'); title(['Output (' label ')']); xlim([0 19.99]); grid on;
    subplot(3,4,row_start+2);
    plot(t, x_real(:,1), 'b-', t, x_hat(:,1), 'r--'); xlabel('Time'); ylabel('x_1'); title(['State 1 (' label ')']); legend('Real','Est'); xlim([0 19.99]); grid on;
    subplot(3,4,row_start+3);
    plot(t, x_real(:,2), 'b-', t, x_hat(:,2), 'r--'); xlabel('Time'); ylabel('x_2'); title(['State 2 (' label ')']); legend('Real','Est'); xlim([0 19.99]); grid on;
end

function [x_real, x_hat, u] = simulate_robustness(steps, dt, t, x_ss, A, B, C, u_nom, K, Lp, Li, x1,x20,x30,x40,x60,c0,c1,c2,c3,c4, c0_new,c1_new,c2_new,c3_new,c4_new, use_saturation)
    u_min = 0; u_max = 2;
    x_real = zeros(steps,2); x_hat = zeros(steps,2); z = zeros(steps,1); u = zeros(steps,1);
    x_real(1,:) = x_ss'; x_hat(1,:) = x_ss';
    for k = 1:steps-1
        x2 = x_real(k,1); x3 = x_real(k,2);
        x_hat_k = x_hat(k,:)';
        d_val = 0.05 * (t(k) >= 5);
        u_hat = u_nom - K(1)*(x_hat_k(1)-x_ss(1)) - K(2)*(x_hat_k(2)-x_ss(2));
        if use_saturation
            if u_hat > u_max, u_hat = u_max; elseif u_hat < u_min, u_hat = u_min; end
        end
        u(k) = u_hat;
        u_real = u_hat + d_val;
        % Real plant with perturbed parameters
        dx2_real = (x60/x1)*(x20-x2) - c0_new*x2*exp(-1/x3);
        dx3_real = (x60/x1)*(x30-x3) + c1_new*x2*exp(-1/x3) - (c2_new*c3_new*u_real*(x3-x40))/(x1*(c3_new*u_real+c4_new));
        y_meas = C * [x2; x3];
        error_y = y_meas - C * x_hat_k;
        xh1 = x_hat_k(1); xh2 = x_hat_k(2);
        % Observer still uses nominal parameters
        dxh1 = (x60/x1)*(x20-xh1) - c0*xh1*exp(-1/xh2);
        dxh2 = (x60/x1)*(x30-xh2) + c1*xh1*exp(-1/xh2) - (c2*c3*u_hat*(xh2-x40))/(x1*(c3*u_hat+c4));
        dx_hat = [dxh1; dxh2] + Lp*error_y + Li*z(k);
        x_real(k+1,:) = x_real(k,:) + dt*[dx2_real,dx3_real];
        x_hat(k+1,:) = x_hat(k,:) + dt*dx_hat';
        z(k+1) = z(k) + dt*error_y;
    end
    u(end) = u(end-1);
end