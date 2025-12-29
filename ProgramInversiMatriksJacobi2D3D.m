clc; clear all; close all;

%% DEFINISI parameter model gravitasi
kons_G = 6.67*10^-11;
R = 70;
rho_model = 100;

% Pusat model bawah permukaan
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

%% Plot 2D Data Observasi dan Model Asli
figure(1)
subplot(2,2,1)
plot(xlin, g_obs, 'or', 'MarkerFaceColor', 'r');
xlabel('Jarak (m)'); ylabel('Gravity Anomaly (mGal)');
title('Grafik g Observasi');

subplot(2,2,3)
rectangle('Position', [x0-R, z0-R, 2*R, 2*R], 'Curvature', [1,1], 'FaceColor', 'r');
daspect([1,1,1]); ylim([0,500]); xlim([0,1000]);
set(gca, 'ydir', 'reverse');
xlabel('Jarak (m)'); ylabel('Kedalaman (m)');
title('Model Bawah Permukaan');

%% 2D - Prediksi Model Awal dan Inversi Jacobi
iterasi = 1;
eps = 1;
e_plot = [];

% Setting plot Model Prediksi
subplot(2,2,4)
hold on
daspect([1,1,1]); ylim([0,500]); xlim([0,1000]);
set(gca, 'ydir', 'reverse');
xlabel('Jarak (m)'); ylabel('Kedalaman (m)');
title('Model Prediksi');

lingkaran_penuh_sudah = false; 

while eps >= 0.0005
    if iterasi == 1
        x0_model = 100;
        z0_model = 100;
    else
        x0_model = x0_pertu;
        z0_model = z0_pertu;
    end
    
    % g_cal prediksi
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
        for i = 1:length(xlin)
            derivative_x(i) = kons_G * (4/3 * pi * R^3 * rho_model) *...
                (3 * z0_model * (xlin(i) - x0_model)) * 10^5 /...
                ((xlin(i)-x0_model)^2 + (zlin(i)-z0_model)^2)^(5/2);
            derivative_z(i) = kons_G * (4/3 * pi * R^3 * rho_model) *...
                (zlin(i)^2 + zlin(i)*z0_model + xlin(i)^2 - 2*xlin(i)*...
                x0_model - 2*z0_model^2 + x0_model^2) * 10^5 / ((xlin(i)-x0_model)^2 + (zlin(i)-z0_model)^2)^(5/2);
        end
        
        J = [derivative_x derivative_z];
        dm_perturbasi = inv(J' * J) * J' * dg_misfit;
        x0_pertu = x0_model + dm_perturbasi(1);
        z0_pertu = z0_model + dm_perturbasi(2);

        jarak = sqrt((x0_pertu - x0)^2 + (z0_pertu - z0)^2);

        if (jarak <= 50) && (~lingkaran_penuh_sudah)
            rectangle('Position', [x0_pertu-R, z0_pertu-R, 2*R, 2*R], ...
                      'Curvature', [1,1], 'FaceColor', 'b', 'EdgeColor', 'b');
            lingkaran_penuh_sudah = true; 
        else
            rectangle('Position', [x0_pertu-R, z0_pertu-R, 2*R, 2*R], ...
                      'Curvature', [1,1], 'EdgeColor', 'b', 'LineWidth', 2);
        end
    end
    iterasi = iterasi + 1;
end

subplot(2,2,2)
plot(xlin, g_cal, 'ob', 'MarkerFaceColor', 'b'); hold on;
plot(xlin, g_obs, '-r');
xlabel('Jarak (m)'); ylabel('Gravity Anomaly (mGal)');
title('Grafik g Kalkulasi');

%% Plot Grafik Misfit 2D
figure(2)
plot(1:length(e_plot), e_plot, 'm-', 'LineWidth', 2)
xlabel('Iterasi'); ylabel('Std Misfit');
title('Grafik Misfit (2D)');
grid on
ylim([0 inf])


%% 3D - Data Observasi
g_obs_3D = zeros(size(X));
for i = 1:length(xlin)
    for j = 1:length(ylin)
        g_obs_3D(j,i) = kons_G * (4/3 * pi * R^3 * z0 * rho_model) * 10^5 / ((X(j,i)-x0)^2 + (Y(j,i)-y0)^2 + z0^2)^(3/2);
    end
end

%% Plot 3D Data Observasi
figure(3);
surf(X, Y, g_obs_3D, 'EdgeColor', 'none');
colormap('jet');
colorbar;
xlabel('Jarak X (m)');
ylabel('Jarak Y (m)');
zlabel('Gravity Anomaly (mGal)');
title('Grafik g Observasi dalam 3D');

%% Plot 3D Model Asli
figure(4);
hold on
axis equal
grid on

x_mesh = linspace(0, 950, 19); 
y_mesh = linspace(0, 950, 19);
z_mesh = linspace(0, 500, 10);

[Xg, Yg, Zg] = meshgrid(x_mesh, y_mesh, z_mesh);

% Garis vertikal (X dan Y fix, Z dari bawah ke atas)
for i = 1:length(x_mesh)
    for j = 1:length(y_mesh)
        plot3([x_mesh(i) x_mesh(i)], [y_mesh(j) y_mesh(j)], [z_mesh(1) z_mesh(end)], 'k');
    end
end

% Garis mendatar X-Z (Y fix)
for i = 1:length(x_mesh)
    for k = 1:length(z_mesh)
        plot3([x_mesh(i) x_mesh(i)], [y_mesh(1) y_mesh(end)], [z_mesh(k) z_mesh(k)], 'k');
    end
end

% Garis mendatar Y-Z (X fix)
for j = 1:length(y_mesh)
    for k = 1:length(z_mesh)
        plot3([x_mesh(1) x_mesh(end)], [y_mesh(j) y_mesh(j)], [z_mesh(k) z_mesh(k)], 'k');
    end
end

[X_model, Y_model, Z_model] = sphere(40); 
X_model = X_model * R + x0;
Y_model = Y_model * R + y0;
Z_model = -Z_model * R + z0;
surf(X_model, Y_model, Z_model, 'FaceAlpha', 0.9, 'EdgeColor', 'none', 'FaceColor', 'r');

plot3(X(:), Y(:), zeros(size(X(:))), 'ko', 'MarkerFaceColor', 'g', 'MarkerSize', 4);

[X_ground, Y_ground] = meshgrid(x_mesh, y_mesh);
Z_ground = zeros(size(X_ground));  
surf(X_ground, Y_ground, Z_ground, 'FaceAlpha', 1, 'EdgeColor', 'none', 'FaceColor', [0.5 0.25 0]); % Coklat

% Label dan pengaturan sumbu
xlabel('Jarak X (m)');
ylabel('Jarak Y (m)');
zlabel('Kedalaman (m)');
title('Model Bawah Permukaan 3D');
set(gca, 'ZDir', 'reverse');

view(45,30);
xlim([0 950]);
ylim([0 950]);
zlim([0 500]);
hold off

%% 3D - Prediksi Model Awal dan Inversi Jacobi
x0_model = 100;
y0_model = 100;
z0_model = 100;
eps = 1;
iterasi = 1;
g_cal_3D = zeros(size(X));
e_plot = [];

while eps >= 0.0005
    for i = 1:length(xlin)
        for j = 1:length(ylin)
            g_cal_3D(j,i) = kons_G * (4/3 * pi * R^3 * z0_model * rho_model) * 10^5 / ((X(j,i)-x0_model)^2 + (Y(j,i)-y0_model)^2 + z0_model^2)^(3/2);
        end
    end

    dg_misfit = g_obs_3D - g_cal_3D;
    eps = std(abs(dg_misfit(:)));
    e_plot(iterasi) = eps;

    derivative_x = zeros(size(X));
    derivative_y = zeros(size(Y));
    derivative_z = zeros(size(X));
    for i = 1:length(xlin)
        for j = 1:length(ylin)
            derivative_x(j,i) = kons_G * (4/3 * pi * R^3 * rho_model) * (3 * z0_model * (X(j,i)-x0_model)) * 10^5 / ((X(j,i)-x0_model)^2 + (Y(j,i)-y0_model)^2 + z0_model^2)^(5/2);
            derivative_y(j,i) = kons_G * (4/3 * pi * R^3 * rho_model) * (3 * z0_model * (Y(j,i)-y0_model)) * 10^5 /((X(j,i)-x0_model)^2 + (Y(j,i)-y0_model)^2 + z0_model^2)^(5/2);
            derivative_z(j,i) = kons_G * (4/3 * pi * R^3 * rho_model) * (2*z0_model - 3*((X(j,i)-x0_model)^2 + (Y(j,i)-y0_model)^2)) * 10^5 /  ((X(j,i)-x0_model)^2 + (Y(j,i)-y0_model)^2 + z0_model^2)^(5/2);
        end
    end
    
    J = [derivative_x(:) derivative_y(:) derivative_z(:)];
    dm_perturbasi = inv(J'*J + 0.01*eye(3)) * J' * dg_misfit(:);
    x0_model = x0_model + dm_perturbasi(1);
    y0_model = y0_model + dm_perturbasi(2);
    z0_model = z0_model + dm_perturbasi(3);

    iterasi = iterasi + 1;
end

%% Plot 3D Model Hasil Inversi
figure(5);
hold on
axis equal
grid on

% Mesh tanah
surf(X_ground, Y_ground, Z_ground, 'FaceAlpha', 1, 'EdgeColor', 'none', 'FaceColor', [0.5 0.25 0]); % Coklat

% Grid garis
for i = 1:length(x_mesh)
    for j = 1:length(y_mesh)
        plot3([x_mesh(i) x_mesh(i)], [y_mesh(j) y_mesh(j)], [z_mesh(1) z_mesh(end)], 'k');
    end
end
for i = 1:length(x_mesh)
    for k = 1:length(z_mesh)
        plot3([x_mesh(i) x_mesh(i)], [y_mesh(1) y_mesh(end)], [z_mesh(k) z_mesh(k)], 'k');
    end
end
for j = 1:length(y_mesh)
    for k = 1:length(z_mesh)
        plot3([x_mesh(1) x_mesh(end)], [y_mesh(j) y_mesh(j)], [z_mesh(k) z_mesh(k)], 'k');
    end
end

% Gambar bola model asli (Merah)
[X_model, Y_model, Z_model] = sphere(40);
X_model = X_model * R + x0;
Y_model = Y_model * R + y0;
Z_model = -Z_model * R + z0;
surf(X_model, Y_model, Z_model, 'FaceAlpha', 0.9, 'EdgeColor', 'none', 'FaceColor', 'r');

% Gambar bola hasil inversi (Biru)
[X_model_inv, Y_model_inv, Z_model_inv] = sphere(40);
X_model_inv = X_model_inv * R + x0_model; % x0_model hasil inversi!
Y_model_inv = Y_model_inv * R + y0_model;
Z_model_inv = -Z_model_inv * R + z0_model;
surf(X_model_inv, Y_model_inv, Z_model_inv, 'FaceAlpha', 0.9, 'EdgeColor', 'none', 'FaceColor', 'b');

% Titik pengukuran
plot3(X(:), Y(:), zeros(size(X(:))), 'ko', 'MarkerFaceColor', 'g', 'MarkerSize', 4);

xlabel('Jarak X (m)');
ylabel('Jarak Y (m)');
zlabel('Kedalaman (m)');
title('Model Bawah Permukaan Hasil Inversi 3D');
set(gca, 'ZDir', 'reverse');
view(45,30);
xlim([0 950]);
ylim([0 950]);
zlim([0 500]);
hold off

%% Plot Grafik Misfit 3D
figure(6);
plot(1:length(e_plot), e_plot, '-m', 'LineWidth', 2);
xlabel('Iterasi');
ylabel('Std Misfit');
title('Grafik Misfit (3D)');

%% Plot Perbandingan Observasi dan Kalkulasi 3D
figure(7);
hold on;
surf(X, Y, g_obs_3D, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
surf(X, Y, g_cal_3D, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
colormap('jet');
xlabel('Jarak X (m)');
ylabel('Jarak Y (m)');
zlabel('Gravity Anomaly (mGal)');
title('Perbandingan g Observasi dan g Kalkulasi dalam 3D');
legend({'g Observasi', 'g Kalkulasi'});
