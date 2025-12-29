clc; clear; close all;

%% Konstanta & model tetap
G   = 6.67e-11;
R   = 70;

% Posisi sumber DIANGGAP BENAR
x0 = 500;
z0 = 250;

% Tebakan awal densitas
rho_model = 80;

%% Lintasan pengukuran
xlin = 0:50:950;
zlin = zeros(size(xlin));

%% Data observasi (? asli = 100)
rho_true = 100;
g_obs = zeros(length(xlin),1);
for i = 1:length(xlin)
    g_obs(i) = G*(4/3*pi*R^3*z0*rho_true)*1e5 / ...
        ((xlin(i)-x0)^2 + (zlin(i)-z0)^2)^(3/2);
end

%% Inversi densitas (Jacobi 1-parameter)
eps = 1;
iter = 1;
lambda = 0.01;

while eps > 5e-4
    
    % Forward model
    g_cal = zeros(length(xlin),1);
    for i = 1:length(xlin)
        g_cal(i) = G*(4/3*pi*R^3*z0*rho_model)*1e5 / ...
            ((xlin(i)-x0)^2 + (zlin(i)-z0)^2)^(3/2);
    end
    
    % Misfit
    dg = g_obs - g_cal;
    eps = std(abs(dg));
    
    % Jacobian (?g/??)
    J = zeros(length(xlin),1);
    for i = 1:length(xlin)
        J(i) = G*(4/3*pi*R^3*z0)*1e5 / ...
            ((xlin(i)-x0)^2 + (zlin(i)-z0)^2)^(3/2);
    end
    
    % Update densitas
    drho = (J'*J + lambda)\(J'*dg);
    rho_model = rho_model + drho;
    
    iter = iter + 1;
end

fprintf('Estimated density (rho): %.2f kg/m^3\n', rho_model);
