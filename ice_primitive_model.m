% ice primitive model
clear; close all; clc;

% configuration
Lx = 120000; Lz = 62000;
nx = 150; nz = 60;
t_end = 300000;
target_cfl = 0.95;
save_path = '/Users/victornguyen/Documents/EMES415/ice_flow/parameter sweep/'; 

% ice sheet geometry and physics
ice_width = Lx;
ice_thickness = 31000;

% constants
g = 0.113; alpha = 5e-4; 
nu = 0.5; kappa = 0.005; 
T_ocean = 0.0; T_ice = -10.0;
T_seafloor = 45.0;      

% grid
dx = Lx/nx; dz = Lz/nz;
x_c = linspace(dx/2, Lx-dx/2, nx);
z_c = linspace(dz/2, Lz-dz/2, nz);
[Xc, Zc] = meshgrid(x_c, z_c);

% initialize fields
T = T_ocean * ones(nz, nx) + 0.01*randn(nz, nx);
p = zeros(nz, nx);
u = zeros(nz, nx);
w = zeros(nz+1, nx);

% storage
T_history = zeros(nz, nx, 100);
u_history = zeros(nz, nx, 100);
w_history = zeros(nz+1, nx, 100);
time_history = zeros(100, 1);
save_interval = t_end / 100;
next_save = 0;
save_step = 1;

% pressure poisson matrix
N = nx * nz;
L = sparse(N, N);
for i = 1:nx
    for j = 1:nz
        row = (i-1)*nz + j;
        left = i-1; if left < 1, left = nx; end
        right = i+1; if right > nx, right = 1; end
        L(row, (left-1)*nz + j) = 1/dx^2;
        L(row, (right-1)*nz + j) = 1/dx^2;
        L(row, row) = L(row, row) - 2/dx^2;
        if j > 1,  L(row, row) = L(row, row) - 1/dz^2; L(row, row-1) = 1/dz^2;  end
        if j < nz, L(row, row) = L(row, row) - 1/dz^2; L(row, row+1) = 1/dz^2;  end
        if i == 1 && j == 1, L(row, :) = 0; L(row, row) = 1; end
    end
end

% main solver loop
t = 0;
fprintf('simulating...\n');
while t < t_end
    % calculate adaptive dt
    max_u = max(abs(u(:)));
    max_w = max(abs(w(:)));
    vel_max = max([max_u, max_w, 1e-4]);
    dt = min(target_cfl * min(dx, dz) / vel_max, 2000.0);
    
    % ice boundary
    ice_mask = (Zc <= ice_thickness);
    T(ice_mask) = T_ice; 
    
    % force zero velocity in ice
    u(Zc <= ice_thickness) = 0;
    w(linspace(0, Lz, nz+1)' <= ice_thickness, :) = 0;
    
    % advection-diffusion
    uc = 0.5 * (u + u(:, [2:nx, 1]));
    wc = 0.5 * (w(1:nz, :) + w(2:nz+1, :));
    
    T_L = T(:, [nx, 1:nx-1]); T_R = T(:, [2:nx, 1]);
    T_U = [T(1, :); T(1:nz-1, :)]; T_D = [T(2:nz, :); T(nz, :)];
    
    adv_x = max(uc,0).*(T-T_L)/dx + min(uc,0).*(T_R-T)/dx;
    adv_z = max(wc,0).*(T-T_U)/dz + min(wc,0).*(T_D-T)/dz;
    diff = kappa * ((T_R - 2*T + T_L)/dx^2 + (T_D - 2*T + T_U)/dz^2);
    
    T(2:nz-1, :) = T(2:nz-1, :) + dt * (- (adv_x(2:nz-1, :) + adv_z(2:nz-1, :)) + diff(2:nz-1, :));
    T(nz, :) = T_seafloor;
    
    % intermediate velocities
    u_star = u; w_star = w;
    u_L = u(:, [nx, 1:nx-1]); u_R = u(:, [2:nx, 1]);
    u_U = [u(1, :); u(1:nz-1, :)]; u_D = [u(2:nz, :); u(nz, :)];
    w_at_u = 0.25*(w(1:nz, :) + w(2:nz+1, :) + w(1:nz, [nx, 1:nx-1]) + w(2:nz+1, [nx, 1:nx-1]));
    
    adv_ux = max(u,0).*(u-u_L)/dx + min(u,0).*(u_R-u)/dx;
    adv_uz = max(w_at_u,0).*(u-u_U)/dz + min(w_at_u,0).*(u_D-u)/dz;
    diff_u = nu * ((u_R - 2*u + u_L)/dx^2 + (u_D - 2*u + u_U)/dz^2);
    u_star(2:nz-1, :) = u(2:nz-1, :) + dt * (-(adv_ux(2:nz-1, :) + adv_uz(2:nz-1, :)) + diff_u(2:nz-1, :));
    
    w_L = w(:, [nx, 1:nx-1]); w_R = w(:, [2:nx, 1]);
    w_U = [w(1, :); w(1:nz, :)]; w_D = [w(2:nz+1, :); w(nz+1, :)];
    
    u_pad = [zeros(1, nx); u; zeros(1, nx)];
    u_at_w_Z = 0.5 * (u_pad(1:nz+1, :) + u_pad(2:nz+2, :));
    u_at_w = 0.5 * (u_at_w_Z + u_at_w_Z(:, [2:nx, 1]));
    
    adv_wx = max(u_at_w,0).*(w-w_L)/dx + min(u_at_w,0).*(w_R-w)/dx;
    adv_wz = max(w,0).*(w-w_U)/dz + min(w,0).*(w_D-w)/dz;
    diff_w = nu * ((w_R - 2*w + w_L)/dx^2 + (w_D - 2*w + w_U)/dz^2);
    buoyancy = g * alpha * (T_ocean - 0.5*(T([1, 1:nz-1], :) + T));
    w_star(2:nz, :) = w(2:nz, :) + dt * (-(adv_wx(2:nz, :) + adv_wz(2:nz, :)) + diff_w(2:nz, :) + buoyancy(2:nz, :));
    
    % pressure projection
    div = (u(:, [2:nx, 1]) - u)/dx + (w(2:nz+1, :) - w(1:nz, :))/dz;
    rhs = div(:) / dt; rhs(1) = 0;
    p_vec = L \ rhs; p = reshape(p_vec, nz, nx);
    
    u = u_star - dt * (p - p(:, [nx, 1:nx-1]))/dx;
    w(2:nz, :) = w_star(2:nz, :) - dt * (p(2:nz, :) - p(1:nz-1, :))/dz;
    t = t + dt;
    
    % progress and save
    if t >= next_save && save_step <= 100
        T_history(:,:,save_step) = T;
        u_history(:,:,save_step) = u;
        w_history(:,:,save_step) = w;
        time_history(save_step) = t;
        save_step = save_step + 1;
        next_save = next_save + save_interval;
    end
end

% visualization and gif export
gif_filename = [save_path, 'ice_diagnostics_sim.gif'];
T_min_global = T_ice;
T_max_global = T_seafloor;
KE_max_global = 20;
Vort_max_global = 0.01;

fig = figure('Color', 'w', 'Position', [50 100 1200 800]);
for k = 1:save_step-1
    if ~ishandle(fig), break; end
    uk = u_history(:,:,k); wk = w_history(:,:,k); Tk = T_history(:,:,k);
    
    subplot(3, 1, 1);
    imagesc(x_c/1000, z_c/1000, Tk); hold on;
    colormap(gca, sky); colorbar; clim([T_min_global, T_max_global]);
    title(sprintf('Enceladus Circulation | t = %.0f s', time_history(k)));
    ylabel('Depth (km)'); set(gca, 'YDir', 'reverse');
    
    subplot(3, 1, 2);
    uc_v = 0.5 * (uk + uk(:, [2:nx, 1]));
    wc_v = 0.5 * (wk(1:nz, :) + wk(2:nz+1, :));
    vorticity = (wc_v(:, [2:nx, 1]) - wc_v(:, [nx, 1:nx-1]))/(2*dx) - ...
                (uc_v([2:nz, nz], :) - uc_v([1, 1:nz-1], :))/(2*dz);
    imagesc(x_c/1000, z_c/1000, vorticity);
    colormap(gca, redblue); colorbar; clim([-Vort_max_global, Vort_max_global]);
    title('Vorticity (1/s)');
    ylabel('Depth (km)'); set(gca, 'YDir', 'reverse');
    
    subplot(3, 1, 3);
    uc_sq = 0.5 * (uk.^2 + uk(:, [2:nx, 1]).^2);
    wc_sq = 0.5 * (wk(1:nz, :).^2 + wk(2:nz+1, :).^2);
    ke = 0.5 * (uc_sq + wc_sq);
    imagesc(x_c/1000, z_c/1000, ke);
    colormap(gca, parula); colorbar; clim([0, KE_max_global]);
    title('Kinetic Energy (m^2/s^2)');
    xlabel('Distance (km)'); ylabel('Depth (km)'); set(gca, 'YDir', 'reverse');
    
    drawnow; frame = getframe(fig); im = frame2im(frame); [imind, cm] = rgb2ind(im, 256);
    if k == 1
        imwrite(imind, cm, gif_filename, 'gif', 'Loopcount', inf, 'DelayTime', 0.1);
    else
        imwrite(imind, cm, gif_filename, 'gif', 'WriteMode', 'append', 'DelayTime', 0.1);
    end
end

% diagnostic plots
fig2 = figure('Color', 'w', 'Position', [100 100 1100 350]);
subplot(1, 3, 1); 
z_c_vec = linspace(dz/2, Lz-dz/2, nz)';
liquid_mask = (z_c_vec > ice_thickness);
T_liquid_mean = zeros(save_step-1, 1);
for k = 1:save_step-1
    Tk = T_history(:,:,k);
    T_liquid_mean(k) = mean(Tk(liquid_mask, :), 'all');
end
plot(time_history(1:save_step-1)/3600, T_liquid_mean, 'b', 'LineWidth', 2);
grid on; xlabel('Time (hours)'); ylabel('Mean Temp (C)'); title('Liquid Ocean Heat');

subplot(1, 3, 2); 
ke_avg = zeros(save_step-1, 1);
for k = 1:save_step-1
    uk = u_history(:,:,k); wk = w_history(:,:,k);
    uc = 0.5 * (uk + uk(:, [2:nx, 1])); 
    wc = 0.5 * (wk(1:nz, :) + wk(2:nz+1, :));
    ke_2d = 0.5 * (uc.^2 + wc.^2);
    ke_avg(k) = mean(ke_2d(liquid_mask, :), 'all');
end
plot(time_history(1:save_step-1)/3600, ke_avg, 'r', 'LineWidth', 2);
grid on; xlabel('Time (hours)'); ylabel('Mean KE'); title('Liquid Ocean KE');

subplot(1, 3, 3);
T_final = T_history(:,:,save_step-1);
T_profile = mean(T_final, 2);
plot(T_profile, z_c/1000, 'k', 'LineWidth', 2);
grid on; xlabel('Temp (C)'); ylabel('Depth (km)');
title('Vertical Temp Profile');
set(gca, 'YDir', 'reverse');
saveas(fig2, [save_path, 'simulation_diagnostics.png']);
