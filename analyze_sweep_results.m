% analyze sweep results
clear; close all; clc;

% global plot settings
set(0, 'DefaultAxesFontSize', 16);
set(0, 'DefaultTextFontSize', 16);
set(0, 'DefaultLineLineWidth', 2);

save_path = '/Users/victornguyen/Documents/EMES415/ice_flow/parameter sweep/';
if ~exist(save_path, 'dir'), mkdir(save_path); end
data_file = fullfile(save_path, 'sweep_results.mat');
if ~exist(data_file, 'file')
    error('sweep_results.mat not found in %s. Please run parameter_sweep_enceladus.m first.', save_path);
end
load(data_file);
nx = 150; nz = 60;
x_km = linspace(0, 120, nx); z_km = linspace(0, 62, nz);

% define sweeps
sweeps = {
    'ice',      [2, 1, 3], {'15km', '31km (Ref)', '45km'}, 'Ice Thickness';
    'seafloor', [4, 5, 1, 6], {'0.0C (None)', '1.5C (Low)', '45.0C (Ref)', '93.0C (Extreme)'}, 'Seafloor Temp'
};

% define variables for 2d maps
variables = {
    'T_mean_2d',    'Temperature',       sky,      [0, 2.0],   '(C)';
    'Vort_mean_2d', 'Vorticity',         redblue,  [-0.01, 0.01], '(s^{-1})';
    'KE_mean_2d',   'Kinetic Energy',   parula,   [0, 20],  '(m^2/s^2)'
};

% 2d comparison figures
for s = 1:size(sweeps, 1)
    sweep_name = sweeps{s, 1}; case_indices = sweeps{s, 2}; case_labels = sweeps{s, 3}; sweep_title = sweeps{s, 4};
    for v = 1:size(variables, 1)
        var_field = variables{v, 1}; var_cmap = variables{v, 3}; var_clim = variables{v, 4}; var_unit = variables{v, 5};
        fig_name = sprintf('fig_%s_sweep_%s', sweep_name, var_field(1:4));
        num_cols = length(case_indices);
        figure('Visible', 'off', 'Color', 'w', 'Position', [50 50 800*num_cols 800]);
        for p = 1:num_cols
            idx = case_indices(p); res = all_results.(sprintf('case_%d', idx));
            subplot(1, num_cols, p); imagesc(x_km, z_km, res.(var_field));
            
            % adjust color limits
            current_clim = var_clim;
            if strcmp(var_field, 'T_mean_2d')
                ts_case = res.T_seafloor; 
                if ts_case <= 0, ts_case = 1.0; end % safety for no-heat case
                current_clim = [-10, 45];
            elseif strcmp(var_field, 'KE_mean_2d')
                current_clim = [0, 50];
            elseif strcmp(var_field, 'Vort_mean_2d')
                current_clim = [-0.01, 0.01];
            end
            
            colormap(gca, var_cmap); cb = colorbar; clim(current_clim);
            ylabel(cb, var_unit);
            title(case_labels{p}); ylabel('Depth (km)'); xlabel('Dist (km)');
            set(gca, 'YDir', 'reverse'); ylim([0, 62]);
        end
        saveas(gcf, fullfile(save_path, [fig_name, '.png'])); close(gcf);
    end
end

% quantitative profiles
for s = 1:size(sweeps, 1)
    sweep_name = sweeps{s, 1}; case_indices = sweeps{s, 2}; case_labels = sweeps{s, 3}; sweep_title = sweeps{s, 4};
    
    % thermal profile
    figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1000 800]); hold on;
    other_colors = {'b', 'r', 'g', 'm'}; color_idx = 1;
    for p = 1:length(case_indices)
        res = all_results.(sprintf('case_%d', case_indices(p)));
        label = case_labels{p};
        if contains(label, '(Ref)')
            line_style = 'k-'; lw = 4;
        else
            line_style = [other_colors{mod(color_idx-1, length(other_colors))+1}, '--'];
            lw = 2; color_idx = color_idx + 1;
        end
        plot(res.T_profile, z_km, line_style, 'LineWidth', lw); 
    end
    grid on; set(gca, 'YDir', 'reverse'); ylim([0, 62]);
    title(sprintf('Temperature vs %s', sweep_title)); xlabel('Temperature (C)'); ylabel('Depth (km)');
    legend(case_labels, 'Location', 'best');
    saveas(gcf, fullfile(save_path, sprintf('quant_thermal_%s.png', sweep_name))); close(gcf);
    
    % kinetic energy profile
    figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1000 800]); hold on;
    color_idx = 1;
    for p = 1:length(case_indices)
        res = all_results.(sprintf('case_%d', case_indices(p)));
        ke_prof = mean(res.KE_mean_2d, 2);
        label = case_labels{p};
        if contains(label, '(Ref)')
            line_style = 'k-'; lw = 4;
        else
            line_style = [other_colors{mod(color_idx-1, length(other_colors))+1}, '--'];
            lw = 2; color_idx = color_idx + 1;
        end
        plot(ke_prof, z_km, line_style, 'LineWidth', lw); 
    end
    grid on; set(gca, 'YDir', 'reverse'); ylim([0, 62]);
    title(sprintf('Kinetic Energy vs %s', sweep_title)); xlabel('Kinetic Energy (m^2/s^2)'); ylabel('Depth (km)');
    legend(case_labels, 'Location', 'best');
    saveas(gcf, fullfile(save_path, sprintf('quant_energy_%s.png', sweep_name))); close(gcf);
end

fprintf('analysis complete.\n');
