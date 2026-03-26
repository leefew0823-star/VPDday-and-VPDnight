clear; clc; close all;

filename = '---'; 
data = readmatrix(filename,'Sheet','Sheet1'); 

num_sites = 11; 

universal_vpd = linspace(0, 1, 100)'; 
export_data = zeros(100, num_sites + 1);
export_data(:, 1) = universal_vpd;

summary_table = cell(num_sites + 1, 13);
summary_table(1, :) = {'Site', 'Best_Model', 'Delta_AIC', 'Adj_R2', 'RMSE', ...
                       'Model_P_val', 'Breakpoint_VPD', ...
                       'Slope_1 (k1)', 'P_val_k1', ...
                       'Slope_2 (k2)', 'P_val_k2', ...
                       'Linear_Slope (k)', 'P_val_k'};

figure('Position', [100, 100, 950, 600]); 
hold on; grid on;

colors = parula(num_sites); 
h_lines = zeros(1, num_sites);
legend_labels = cell(1, num_sites);

for i = 1:num_sites
    col_vpd = 2 * i - 1; 
    col_sap = 2 * i;     
    
    vpd = data(:, col_vpd);
    sap = data(:, col_sap);
    
    valid_idx = ~isnan(vpd) & ~isnan(sap);
    vpd_valid = vpd(valid_idx);
    sap_valid = sap(valid_idx);
    
   
    vpd_min = min(vpd_valid); vpd_max = max(vpd_valid);
    if vpd_max > vpd_min; vpd_valid = (vpd_valid - vpd_min) / (vpd_max - vpd_min);
    else; vpd_valid = zeros(size(vpd_valid)); end

    sap_min = min(sap_valid); sap_max = max(sap_valid);
    if sap_max > sap_min; sap_valid = (sap_valid - sap_min) / (sap_max - sap_min);
    else; sap_valid = zeros(size(sap_valid)); end
    
    scatter(vpd_valid, sap_valid, 15, colors(i,:), 'o', 'filled', 'MarkerFaceAlpha', 0.15);
    
    [x_sort, sort_idx] = sort(vpd_valid);
    y_sort = sap_valid(sort_idx);
    n_points = length(x_sort);
    
    SST = sum((y_sort - mean(y_sort)).^2);
    
    % Linear
    p_lin = polyfit(x_sort, y_sort, 1);
    y_fit_lin = polyval(p_lin, x_sort);
    SSR_lin = sum((y_sort - y_fit_lin).^2);
    
    AIC_lin = n_points * log(SSR_lin / n_points + 1e-10) + 2 * 2;
    R2_lin = 1 - (SSR_lin / SST);
    Adj_R2_lin = 1 - ((1 - R2_lin) * (n_points - 1)) / (n_points - 2); 
    RMSE_lin = sqrt(SSR_lin / n_points);

    F_lin = ((SST - SSR_lin) / 1) / (SSR_lin / (n_points - 2));
    p_val_model_lin = 1 - fcdf(F_lin, 1, n_points - 2);
   
    SSx_lin = sum((x_sort - mean(x_sort)).^2);
    SE_k = sqrt((SSR_lin / (n_points - 2)) / SSx_lin); 
    t_k = p_lin(1) / SE_k; 
    p_val_k = 1 - fcdf(t_k^2, 1, n_points - 2); 
    
    % Piecewise
    best_SSR_pw = inf;      
    best_bp = x_sort(1); 
    best_p1 = [0, 0];    
    best_k2 = 0;         
    best_idx = 1;
    
    min_idx = max(3, round(n_points * 0.1));
    max_idx = min(n_points - 3, round(n_points * 0.9));
    
    for j = min_idx:max_idx
        bp = x_sort(j); 
        x1 = x_sort(1:j);     y1 = y_sort(1:j);
        x2 = x_sort(j+1:end); y2 = y_sort(j+1:end);
        
        p1 = polyfit(x1, y1, 1);
        y1_fit = polyval(p1, x1);
        y_at_bp = polyval(p1, bp);
        
        if ~isempty(x2)
            k2 = (x2 - bp) \ (y2 - y_at_bp); 
            y2_fit = k2 * (x2 - bp) + y_at_bp;
        else
            k2 = 0;
            y2_fit = [];
        end
        
        SSR_pw = sum((y1 - y1_fit).^2) + sum((y2 - y2_fit).^2);
        
        if SSR_pw < best_SSR_pw
            best_SSR_pw = SSR_pw;
            best_bp = bp;
            best_p1 = p1;
            best_k2 = k2;
            best_idx = j; 
        end
    end
    
    AIC_pw = n_points * log(best_SSR_pw / n_points + 1e-10) + 2 * 4;
    R2_pw = 1 - (best_SSR_pw / SST);
    Adj_R2_pw = 1 - ((1 - R2_pw) * (n_points - 1)) / (n_points - 4); 
    RMSE_pw = sqrt(best_SSR_pw / n_points);
    
   
    F_pw = ((SST - best_SSR_pw) / 3) / (best_SSR_pw / (n_points - 4));
    p_val_model_pw = 1 - fcdf(F_pw, 3, n_points - 4);

    x1 = x_sort(1:best_idx);     y1 = y_sort(1:best_idx);
    x2 = x_sort(best_idx+1:end); y2 = y_sort(best_idx+1:end);
    n1 = length(x1);             n2 = length(x2);
    
    % k1
    y1_fit = polyval(best_p1, x1);
    SSR1 = sum((y1 - y1_fit).^2);
    SSx1 = sum((x1 - mean(x1)).^2);
    SE_k1 = sqrt((SSR1 / (n1 - 2)) / SSx1);
    t_k1 = best_p1(1) / SE_k1;
    p_val_k1 = 1 - fcdf(t_k1^2, 1, n1 - 2);
    
    % k2 
    y_at_bp = polyval(best_p1, best_bp);
    y2_fit = best_k2 * (x2 - best_bp) + y_at_bp;
    SSR2 = sum((y2 - y2_fit).^2);
    SSx2 = sum((x2 - best_bp).^2);
    SE_k2 = sqrt((SSR2 / (n2 - 1)) / SSx2);
    t_k2 = best_k2 / SE_k2;
    p_val_k2 = 1 - fcdf(t_k2^2, 1, n2 - 1);
    
    Delta_AIC = AIC_lin - AIC_pw;
 
    y_universal = zeros(100, 1);
    summary_table{i+1, 1} = ['Site ', num2str(i)];
    summary_table{i+1, 3} = round(Delta_AIC, 2);
    
    if Delta_AIC > 2
        summary_table{i+1, 2} = 'Piecewise';
        summary_table{i+1, 4} = round(Adj_R2_pw, 3);
        summary_table{i+1, 5} = round(RMSE_pw, 4);
        summary_table{i+1, 6} = format_p(p_val_model_pw);
        summary_table{i+1, 7} = round(best_bp, 3);
        
        summary_table{i+1, 8} = round(best_p1(1), 3); % k1
        summary_table{i+1, 9} = format_p(p_val_k1);   % p_k1
        
        summary_table{i+1, 10} = round(best_k2, 3);   % k2
        summary_table{i+1, 11} = format_p(p_val_k2);  % p_k2
        
        summary_table{i+1, 12} = '-';                 
        summary_table{i+1, 13} = '-';                 
        
        idx_before = universal_vpd <= best_bp;
        idx_after = universal_vpd > best_bp;
        y_universal(idx_before) = polyval(best_p1, universal_vpd(idx_before));
        y_universal(idx_after) = best_k2 * (universal_vpd(idx_after) - best_bp) + polyval(best_p1, best_bp);
    else
        summary_table{i+1, 2} = 'Linear';
        summary_table{i+1, 4} = round(Adj_R2_lin, 3);
        summary_table{i+1, 5} = round(RMSE_lin, 4);
        summary_table{i+1, 6} = format_p(p_val_model_lin);
        summary_table{i+1, 7} = '-';
        
        summary_table{i+1, 8} = '-';
        summary_table{i+1, 9} = '-';
        summary_table{i+1, 10} = '-';
        summary_table{i+1, 11} = '-';
        
        summary_table{i+1, 12} = round(p_lin(1), 3);  % k 
        summary_table{i+1, 13} = format_p(p_val_k);   % p_k
        
        y_universal = polyval(p_lin, universal_vpd);
    end
    
    export_data(:, i + 1) = y_universal;
    h_lines(i) = plot(universal_vpd, y_universal, '-', 'Color', colors(i,:), 'LineWidth', 2.5);
    legend_labels{i} = ['Site ', num2str(i)];
end

summary_filename = '---';
writecell(summary_table, summary_filename);

xlabel('Normalized VPD', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Normalized Sap Flow', 'FontSize', 14, 'FontWeight', 'bold'); 
title('Robust Model Fit (Selected via \Delta AIC > 2)', 'FontSize', 14);

legend(h_lines, legend_labels, 'Location', 'eastoutside', 'FontSize', 11);
set(gca, 'FontSize', 12, 'LineWidth', 1.2);
box on;

