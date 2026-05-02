% enceladus_box_convection.m
% 2D Advection-Diffusion: Shrinking Floating Ice on Enceladus
% -------------------------------------------------------------------------

clear; close all; clc;

% 1. Grid & Physical Parameters
Lx = 1000;           % Domain width (m)
Lz = 600;            % Total depth (m)
nx = 151;            
nz = 101;             

dx = Lx / (nx - 1);
dz = Lz / (nz - 1);
[X, Z] = meshgrid(linspace(0, Lx, nx), linspace(0, Lz, nz));

% Physics
g = 0.113;           % Surface gravity (m/s²)
alpha = 2e-4;        % Thermal expansion
nu = 0.05;           % Viscosity
kappa = 0.02;        % Thermal diffusivity
rho_ice = 917;       % kg/m³
rho_water = 1000;    % kg/m³

% Interface & Temperature
z_surface = 150;     % Waterline at 150m
T_bottom = 25.0;     
T_air = 15.0;        
T_ice = -5.0;        
melt_rate = 0.005;   % Shrinking speed factor

% Initial Ice Dimensions
bw = 300;            % Initial width (m)
bd = 150;            % Initial total thickness (m)

% 2. Initial Fields
T = ones(nz, nx) * 5.0;  
T(Z <= z_surface) = T_air;

% Initial Ice Placement (Based on Archimedes)
submerged_depth = bd * (rho_ice / rho_water);
block_x = [Lx/2 - bw/2, Lx/2 + bw/2];
block_z = [z_surface - (bd - submerged_depth), z_surface + submerged_depth];
is_ice = (X >= block_x(1) & X <= block_x(2) & Z >= block_z(1) & Z <= block_z(2));
T(is_ice) = T_ice;

% Dynamics Fields
zeta = zeros(nz, nx);
psi = zeros(nz, nx);
u = zeros(nz, nx); w = zeros(nz, nx);

% Numerical Setup
dt = 2.0;            
t_end = 25000;       
nt = ceil(t_end / dt);
save_interval = 50;  

% Build Poisson Matrix
N = nx * nz;
L_mat = sparse(N, N);
for k = 1:nz
    for j = 1:nx
        idx = (j-1)*nz + k;
        if j > 1 && j < nx && k > 1 && k < nz
            L_mat(idx, idx) = -2/dx^2 - 2/dz^2;
            L_mat(idx, idx-1) = 1/dz^2; L_mat(idx, idx+1) = 1/dz^2;
            L_mat(idx, idx-nz) = 1/dx^2; L_mat(idx, idx+nz) = 1/dx^2;
        else
            L_mat(idx, idx) = 1; 
        end
    end
end

% 3. Main Solver Loop
fprintf('Simulating Shrinking Floating Ice...\n');
fig = figure('Color', 'w', 'Position', [100, 100, 1000, 600]);
gif_filename = '/Users/victornguyen/Documents/EMES415/ice_flow/enceladus_convection_sim.gif';

for n = 1:nt
    % A. SHRINKING LOGIC
    % Ice shrinks proportional to temperature difference in water
    T_ambient = 10.0; % Approximate ambient water temp
    shrink_step = melt_rate * (T_ambient - T_ice) * dt / 1000;
    bw = max(20, bw - shrink_step * 2);
    bd = max(10, bd - shrink_step);
    
    % Update Archimedes floating position
    submerged_depth = bd * (rho_ice / rho_water);
    block_x = [Lx/2 - bw/2, Lx/2 + bw/2];
    block_z = [z_surface - (bd - submerged_depth), z_surface + submerged_depth];
    is_ice = (X >= block_x(1) & X <= block_x(2) & Z >= block_z(1) & Z <= block_z(2));

    % B. ADVECTION-DIFFUSION
    T_old = T;
    T_int = T(2:nz-1, 2:nx-1);
    u_int = u(2:nz-1, 2:nx-1); w_int = w(2:nz-1, 2:nx-1);
    
    adv_T = max(u_int, 0) .* (T_int - T(2:nz-1, 1:nx-2))/dx + ...
            min(u_int, 0) .* (T(2:nz-1, 3:nx) - T_int)/dx + ...
            max(w_int, 0) .* (T_int - T(1:nz-2, 2:nx-1))/dz + ...
            min(w_int, 0) .* (T(3:nz, 2:nx-1) - T_int)/dz;
            
    diff_T = kappa * ((T(2:nz-1, 3:nx) - 2*T_int + T(2:nz-1, 1:nx-2))/dx^2 + ...
                      (T(3:nz, 2:nx-1) - 2*T_int + T(1:nz-2, 2:nx-1))/dz^2);
                      
    T(2:nz-1, 2:nx-1) = T_int + dt * (-adv_T + diff_T);
    
    % Physics Constraints:
    % 1. Air remains undisturbed (Atmosphere is a heat bath)
    T(Z <= z_surface & ~is_ice) = T_air;
    % 2. Seafloor remains heated
    T(end, :) = T_bottom;
    % 3. Ice Block remains cold
    T(is_ice) = T_ice;
    
    % C. DYNAMICS: Buoyancy in Ocean
    is_water = (Z(2:nz-1, 2:nx-1) > z_surface);
    dTdx = (T(2:nz-1, 3:nx) - T(2:nz-1, 1:nx-2)) / (2*dx);
    buoy_torque = g * alpha * dTdx .* is_water;
    
    zeta_old = zeta;
    z_int = zeta(2:nz-1, 2:nx-1);
    adv_z = max(u_int, 0) .* (z_int - zeta(2:nz-1, 1:nx-2))/dx + ...
            min(u_int, 0) .* (zeta(2:nz-1, 3:nx) - z_int)/dx + ...
            max(w_int, 0) .* (z_int - zeta(1:nz-2, 2:nx-1))/dz + ...
            min(w_int, 0) .* (zeta(3:nz, 2:nx-1) - z_int)/dz;
    diff_z = nu * ((zeta(2:nz-1, 3:nx) - 2*z_int + zeta(2:nz-1, 1:nx-2))/dx^2 + ...
                   (zeta(3:nz, 2:nx-1) - 2*z_int + zeta(1:nz-2, 2:nx-1))/dz^2);
                   
    zeta(2:nz-1, 2:nx-1) = z_int + dt * (-adv_z + diff_z + buoy_torque);
    
    psi_vec = L_mat \ (-zeta(:));
    psi = reshape(psi_vec, nz, nx);
    u(2:nz-1, :) = (psi(3:nz, :) - psi(1:nz-2, :)) / (2*dz);
    w(:, 2:nx-1) = -(psi(:, 3:nx) - psi(:, 1:nx-2)) / (2*dx);
    
    % D. Visualization
    if mod(n, save_interval) == 0
        clf;
        imagesc(linspace(0, Lx, nx), linspace(0, Lz, nz), T);
        hold on;
        plot([0 Lx], [z_surface z_surface], 'w--', 'LineWidth', 1.5); % Waterline
        
        if any(is_ice(:))
            plot([block_x(1) block_x(2) block_x(2) block_x(1) block_x(1)], ...
                 [block_z(1) block_z(1) block_z(2) block_z(2) block_z(1)], 'w-', 'LineWidth', 2);
        end
        
        colormap(cmocean('ice')); colorbar; clim([-5 25]);
        title(sprintf('Enceladus: Shrinking Ice & Convection (t = %.1fs)', n*dt));
        xlabel('Distance (m)'); ylabel('Depth (m)');
        set(gca, 'YDir', 'reverse', 'FontSize', 12);
        drawnow;
        
        frame = getframe(fig); im = frame2im(frame); [imind, cm] = rgb2ind(im, 256);
        if n == save_interval
            imwrite(imind, cm, gif_filename, 'gif', 'Loopcount', inf, 'DelayTime', 0.05);
        else
            imwrite(imind, cm, gif_filename, 'gif', 'WriteMode', 'append', 'DelayTime', 0.05);
        end
    end
end
fprintf('Done. Result: %s\n', gif_filename);
