% parameter sweep enceladus
clear; close all; clc;

% define sensitivity experiments
H_ref = 31000; To_ref = 0.0; Ts_ref = 45.0;
save_path = '/Users/victornguyen/Documents/EMES415/ice_flow/parameter sweep/';
if ~exist(save_path, 'dir'), mkdir(save_path); end

% experiment list
exp_list = {
    H_ref, To_ref, Ts_ref, 'reference';
    15000, To_ref, Ts_ref, 'thin_ice';
    45000, To_ref, Ts_ref, 'thick_ice';
    H_ref, To_ref, 0.0,    'no_hydro';
    H_ref, To_ref, 1.5,    'low_hydro';
    H_ref, To_ref, 90.0,   'extreme_hydro'
};

all_results = struct();
num_exps = size(exp_list, 1);
fprintf('starting sweep...\n');

for e = 1:num_exps
    h = exp_list{e, 1}; to = exp_list{e, 2}; ts = exp_list{e, 3}; label = exp_list{e, 4};
    gif_name = fullfile(save_path, sprintf('case_%d_%s_diagnostics.gif', e, label));
    
    % run solver
    [T_prof, KE_final, T_mean_2d, Vort_mean_2d, KE_mean_2d] = run_ice_solver_with_gif(h, to, ts, 300000, gif_name, label);
    
    % store results
    res_name = sprintf('case_%d', e);
    all_results.(res_name).label = label;
    all_results.(res_name).T_profile = T_prof;
    all_results.(res_name).KE = KE_final;
    all_results.(res_name).T_mean_2d = T_mean_2d;
    all_results.(res_name).Vort_mean_2d = Vort_mean_2d;
    all_results.(res_name).KE_mean_2d = KE_mean_2d;
    all_results.(res_name).H_ice = h;
    all_results.(res_name).T_ocean = to;
    all_results.(res_name).T_seafloor = ts;
end

% save results
save(fullfile(save_path, 'sweep_results.mat'), 'all_results');
fprintf('sweep complete.\n');

% helper function: solver + animation
function [T_profile, ke_final, T_avg, Vort_avg, KE_avg] = run_ice_solver_with_gif(H_ice, T_ocean, T_seafloor, t_end, gif_filename, label_text)
    Lx = 120000; Lz = 62000; nx = 150; nz = 60; 
    dx = Lx/nx; dz = Lz/nz;
    g = 0.113; alpha = 5e-4; nu = 0.5; kappa = 0.005; 
    T_ice = -10; target_cfl = 0.95;
    
    T = T_ocean * ones(nz, nx) + 0.01*randn(nz, nx);
    u = zeros(nz, nx); w = zeros(nz+1, nx);
    
    % storage for 50 frames
    num_frames = 50;
    T_hist = zeros(nz, nx, num_frames);
    u_hist = zeros(nz, nx, num_frames);
    w_hist = zeros(nz+1, nx, num_frames);
    t_hist = zeros(num_frames, 1);
    save_interval = t_end / num_frames; next_save = 0; save_idx = 1;
    
    % poisson matrix
    N = nx * nz; L_mat = sparse(N, N);
    for i = 1:nx
        for j = 1:nz
            row = (i-1)*nz + j;
            left = mod(i-2, nx) + 1; right = mod(i, nx) + 1;
            L_mat(row, (left-1)*nz + j) = 1/dx^2;
            L_mat(row, (right-1)*nz + j) = 1/dx^2;
            L_mat(row, row) = L_mat(row, row) - 2/dx^2;
            if j > 1,  L_mat(row, row) = L_mat(row, row) - 1/dz^2; L_mat(row, row-1) = 1/dz^2;  end
            if j < nz, L_mat(row, row) = L_mat(row, row) - 1/dz^2; L_mat(row, row+1) = 1/dz^2;  end
            if i == 1 && j == 1, L_mat(row, :) = 0; L_mat(row, row) = 1; end
        end
    end
    
    z_c = linspace(dz/2, Lz-dz/2, nz)';
    t = 0;
    while t < t_end
        vel_max = max([abs(u(:)); abs(w(:)); 1e-4]);
        dt = min(target_cfl * min(dx, dz) / vel_max, 2000);
        
        T(z_c <= H_ice, :) = T_ice; T(nz, :) = T_seafloor; 
        
        % force zero velocity in ice
        u(z_c <= H_ice, :) = 0;
        w(linspace(0, Lz, nz+1)' <= H_ice, :) = 0;
        
        uc = 0.5 * (u + u(:, [2:nx, 1])); wc = 0.5 * (w(1:nz, :) + w(2:nz+1, :));
        T_L = T(:, [nx, 1:nx-1]); T_R = T(:, [2:nx, 1]);
        T_U = [T(1, :); T(1:nz-1, :)]; T_D = [T(2:nz, :); T(nz, :)];
        adv_T = max(uc,0).*(T-T_L)/dx + min(uc,0).*(T_R-T)/dx + max(wc,0).*(T-T_U)/dz + min(wc,0).*(T_D-T)/dz;
        diff_T = kappa * ((T_R - 2*T + T_L)/dx^2 + (T_D - 2*T + T_U)/dz^2);
        T(2:nz-1, :) = T(2:nz-1, :) + dt * (-adv_T(2:nz-1, :) + diff_T(2:nz-1, :));
        
        % intermediate velocities
        u_star = u; w_star = w;
        
        % u-velocity
        u_L = u(:, [nx, 1:nx-1]); u_R = u(:, [2:nx, 1]);
        u_U = [u(1, :); u(1:nz-1, :)]; u_D = [u(2:nz, :); u(nz, :)];
        w_at_u = 0.25*(w(1:nz, :) + w(2:nz+1, :) + w(1:nz, [nx, 1:nx-1]) + w(2:nz+1, [nx, 1:nx-1]));
        
        adv_ux = max(u,0).*(u-u_L)/dx + min(u,0).*(u_R-u)/dx;
        adv_uz = max(w_at_u,0).*(u-u_U)/dz + min(w_at_u,0).*(u_D-u)/dz;
        diff_u = nu * ((u_R - 2*u + u_L)/dx^2 + (u_D - 2*u + u_U)/dz^2);
        u_star(2:nz-1, :) = u(2:nz-1, :) + dt * (-(adv_ux(2:nz-1, :) + adv_uz(2:nz-1, :)) + diff_u(2:nz-1, :));
        
        % w-velocity
        w_L = w(:, [nx, 1:nx-1]); w_R = w(:, [2:nx, 1]);
        w_U = [w(1, :); w(1:nz, :)]; w_D = [w(2:nz+1, :); w(nz+1, :)];
        
        u_pad = [zeros(1, nx); u; zeros(1, nx)];
        u_at_w_Z = 0.5 * (u_pad(1:nz+1, :) + u_pad(2:nz+2, :));
        u_at_w = 0.5 * (u_at_w_Z + u_at_w_Z(:, [2:nx, 1]));
        
        adv_wx = max(u_at_w,0).*(w-w_L)/dx + min(u_at_w,0).*(w_R-w)/dx;
        adv_wz = max(w,0).*(w-w_U)/dz + min(w,0).*(w_D-w)/dz;
        diff_w = nu * ((w_R - 2*w + w_L)/dx^2 + (w_D - 2*w + w_U)/dz^2);
        buoy = g * alpha * (T_ocean - 0.5*(T(1:nz-1, :) + T(2:nz, :)));
        w_star(2:nz, :) = w(2:nz, :) + dt * (-(adv_wx(2:nz, :) + adv_wz(2:nz, :)) + diff_w(2:nz, :) + buoy);
        
        div = (u_star(:, [2:nx, 1]) - u_star)/dx + (w_star(2:nz+1, :) - w_star(1:nz, :))/dz;
        rhs = div(:) / dt; rhs(1) = 0;
        p = reshape(L_mat \ rhs, nz, nx);
        u = u_star - dt * (p - p(:, [nx, 1:nx-1]))/dx;
        w(2:nz, :) = w_star(2:nz, :) - dt * (p(2:nz, :) - p(1:nz-1, :))/dz;
        
        if t >= next_save && save_idx <= num_frames
            T_hist(:,:,save_idx) = T; u_hist(:,:,save_idx) = u; w_hist(:,:,save_idx) = w;
            t_hist(save_idx) = t; next_save = next_save + save_interval; save_idx = save_idx + 1;
            fprintf('.');
        end
        t = t + dt;
    end
    
    % generate gif
    T_min_g = T_ice; T_max_g = T_seafloor; 
    V_max_g = 0.01; 
    KE_max_g = 20;
    
    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [50 100 800 600]);
    x_km = linspace(0, 120, nx); z_km = linspace(0, 62, nz);
    for k = 1:num_frames
        subplot(3, 1, 1); imagesc(x_km, z_km, T_hist(:,:,k)); colormap(gca, sky); colorbar; clim([T_min_g, T_max_g]);
        clean_label = strrep(label_text, '_', ' ');
        title(sprintf('%s | t = %.0f s', clean_label, t_hist(k))); ylabel('Depth (km)'); set(gca, 'YDir', 'reverse');
        
        subplot(3, 1, 2); 
        uk = u_hist(:,:,k); wk = w_hist(:,:,k);
        uc_v = 0.5 * (uk + uk(:, [2:nx, 1])); wc_v = 0.5 * (wk(1:nz, :) + wk(2:nz+1, :));
        vort = (wc_v(:, [2:nx, 1]) - wc_v(:, [nx, 1:nx-1]))/(2*dx) - (uc_v([2:nz, nz], :) - uc_v([1, 1:nz-1], :))/(2*dz);
        imagesc(x_km, z_km, vort); colormap(gca, redblue); colorbar; clim([-V_max_g, V_max_g]);
        ylabel('Depth'); set(gca, 'YDir', 'reverse'); title('Vorticity');
        
        subplot(3, 1, 3);
        ke = 0.5 * (uc_v.^2 + wc_v.^2);
        imagesc(x_km, z_km, ke); colormap(gca, parula); colorbar; clim([0, KE_max_g]);
        ylabel('Depth'); set(gca, 'YDir', 'reverse'); title('Kinetic Energy'); xlabel('Dist (km)');
        
        drawnow; frame = getframe(fig); im = frame2im(frame); [imind, cm] = rgb2ind(im, 256);
        if k == 1
            imwrite(imind, cm, gif_filename, 'gif', 'Loopcount', inf, 'DelayTime', 0.1);
        else
            imwrite(imind, cm, gif_filename, 'gif', 'WriteMode', 'append', 'DelayTime', 0.1);
        end
    end
    close(fig);
    T_profile = mean(T, 2); 
    ke_final = 0.5 * (mean(u(:).^2) + mean(w(:).^2));
    
    % compute 2d time averages
    steady_idx = floor(num_frames/2):num_frames;
    T_avg = mean(T_hist(:,:,steady_idx), 3);
    liquid_mask = (z_c > H_ice);
    
    vort_hist = zeros(nz, nx, length(steady_idx));
    ke_hist = zeros(nz, nx, length(steady_idx));
    for s = 1:length(steady_idx)
        idx = steady_idx(s);
        uk = u_hist(:,:,idx); wk = w_hist(:,:,idx);
        uc_v = 0.5 * (uk + uk(:, [2:nx, 1])); wc_v = 0.5 * (wk(1:nz, :) + wk(2:nz+1, :));
        vort_hist(:,:,s) = (wc_v(:, [2:nx, 1]) - wc_v(:, [nx, 1:nx-1]))/(2*dx) - (uc_v([2:nz, nz], :) - uc_v([1, 1:nz-1], :))/(2*dz);
        ke_hist(:,:,s) = 0.5 * (uc_v.^2 + wc_v.^2);
    end
    Vort_avg = mean(vort_hist, 3);
    KE_avg = mean(ke_hist, 3);
    
    % final scalar ke
    ke_final = mean(KE_avg(liquid_mask, :), 'all');
    T_profile = mean(T, 2); 
end
