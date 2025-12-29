clc; clear; close all;

%% DEFINISI parameter model gravitasi
kons_G = 6.67*10^-11;
R = 70;
rho_model = 100;

x0 = 500;
y0 = 250;
z0 = 250;

%% Lintasan Pengukuran
xlin = 0:50:950;
ylin = 0:50:950;
[X, Y] = meshgrid(xlin, ylin);
zlin = zeros(size(X));

%% 2D - DATA Observasi
g_obs = zeros(length(xlin), 1);
for i = 1:length(xlin)
    g_obs(i) = kons_G * (4/3 * pi * R^3 * z0 * rho_model) * 10^5 / ((xlin(i)-x0)^2 + (zlin(i)-z0)^2)^(3/2);
end

%% 2D - Prediksi Model Awal dan Inversi Jacobi
iterasi = 1;
eps = 1;
e_plot = [];
while eps >= 0.0005
    if iterasi == 1
        x0_model = 100;
        z0_model = 100;
        rho_model = 80;
    else
        x0_model = x0_pertu;
        z0_model = z0_pertu;
        rho_model = rho_pertu;
    end
    
    g_cal = zeros(length(xlin), 1);
    for i = 1:length(xlin)
        g_cal(i) = kons_G * (4/3 * pi * R^3 * z0_model * rho_model) * 10^5 /...
            ((xlin(i)-x0_model)^2 + (zlin(i)-z0_model)^2)^(3/2);
    end
    
    dg_misfit = g_obs - g_cal;
    eps = std(abs(dg_misfit));
    e_plot(iterasi) = eps;

    if eps >= 0.0005
        derivative_x = zeros(length(xlin),1);
        derivative_z = zeros(length(xlin),1);
        derivative_rho = zeros(length(xlin),1);
        for i = 1:length(xlin)
            derivative_x(i) = kons_G * (4/3 * pi * R^3 * rho_model) *...
                (3 * z0_model * (xlin(i) - x0_model)) * 10^5 / ((xlin(i)-x0_model)^2 + (zlin(i)-z0_model)^2)^(5/2);
            derivative_z(i) = kons_G * (4/3 * pi * R^3 * rho_model) *...
                (zlin(i)^2 + zlin(i)*z0_model + xlin(i)^2 - 2*xlin(i)*x0_model - 2*z0_model^2 + x0_model^2) * 10^5 /...
                ((xlin(i)-x0_model)^2 + (zlin(i)-z0_model)^2)^(5/2);
            derivative_rho(i) = kons_G * (4/3 * pi * R^3 * z0_model) * 10^5 /...
                ((xlin(i)-x0_model)^2 + (zlin(i)-z0_model)^2)^(3/2);
        end
        
        J = [derivative_x derivative_z derivative_rho];
        lambda = 0.01; 
        dm_perturbasi = inv(J' * J + lambda * eye(size(J,2))) * J' * dg_misfit;
        x0_pertu = x0_model + dm_perturbasi(1);
        z0_pertu = z0_model + dm_perturbasi(2);
        rho_pertu = rho_model + dm_perturbasi(3); 
    end
    iterasi = iterasi + 1;
end

fprintf('Estimasi Densitas (rho) akhir: %.2f kg/m^3\n', rho_model);