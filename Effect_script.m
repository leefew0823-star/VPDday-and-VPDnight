%% Ridge Regression Analysis: Daytime vs Nighttime VPD impact on NDVI
clear; clc;
%% step 1
start_year = 1982;
end_year   = 2020;
n_years    = end_year - start_year + 1; 
n_months   = n_years * 12; 

k_param = 5;     
alpha   = 0.01;   
ridge_k = 5;     

doy_to_month_lut = zeros(1, 366);
current_month = 1;
days_accum = 0;
days_per_month = [31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]; % 按闰年算，保证覆盖366
for m = 1:12
    doy_to_month_lut(days_accum+1 : days_accum+days_per_month(m)) = m;
    days_accum = days_accum + days_per_month(m);
end

%% step 2
load('F:\SOS_SG_Dynamic_1982_2022.mat', 'SOS_Matrix');
load('F:\EOS_SG_Dynamic_1982_2022.mat', 'EOS_Matrix');

SOS_Matrix1 = SOS_Matrix(13:291, :, 1:n_years);
EOS_Matrix1 = EOS_Matrix(13:291, :, 1:n_years);

f1 = load('F:\ndvi4g_deseason_1982_2020.mat');  ndvi_all = f1.ndvi4g_deseason_1982_2020; 
f2 = load('\vpd_day_deseason_1982_2020.mat'); vpd_day_all = f2.vpd_day_deseason_1982_2020;
f3 = load('\vpd_night_deseason_1982_2020.mat'); vpd_night_all = f3.vpd_night_deseason_1982_2020;
f4 = load('\t_day_deseason_1982_2020.mat');   t_day_all = f4.t_day_deseason_1982_2020;
f5 = load('\t_night_deseason_1982_2020.mat'); t_night_all = f5.t_night_deseason_1982_2020;
f6 = load('\SMroot_deseason_1982_2020.mat');  swc_all = f6.SMroot_deseason_1982_2020;
f7 = load('\cld_deseason_1982_2020.mat');     cld_all = f7.cld_deseason_1982_2020;

clear f1 f2 f3 f4 f5 f6 f7;

[rows, cols, ~] = size(ndvi_all);

%% step 3
ridge_vpd_ndvi_day   = NaN(rows, cols);
ridge_vpd_ndvi_night = NaN(rows, cols);

%% step 4 Pixel-wise Loop
parfor i = 1:rows
    
    row_day_coef   = NaN(cols, 1);
    row_night_coef = NaN(cols, 1);

    for j = 1:cols
        % Step A: 
        ts_sos = squeeze(SOS_Matrix1(i, j, :)); 
        ts_eos = squeeze(EOS_Matrix1(i, j, :)); 

        if all(isnan(ts_sos))
            continue;
        end

        v_ndvi = squeeze(ndvi_all(i, j, :));
        v_vpdd = squeeze(vpd_day_all(i, j, :));
        v_vpdn = squeeze(vpd_night_all(i, j, :));
        v_td   = squeeze(t_day_all(i, j, :));
        v_tn   = squeeze(t_night_all(i, j, :));
        v_swc  = squeeze(swc_all(i, j, :));
        v_cld  = squeeze(cld_all(i, j, :));

        if all(isnan(v_ndvi))
            continue;
        end

        % Step B
        mask_gs = false(n_months, 1);

        for y = 1:n_years
            sos_d = ts_sos(y);
            eos_d = ts_eos(y);

            if isnan(sos_d) || isnan(eos_d) || sos_d >= eos_d
                continue; 
            end

            m_start = doy_to_month_lut(max(1, min(366, round(sos_d))));
            m_end   = doy_to_month_lut(max(1, min(366, round(eos_d))));

            year_offset = (y - 1) * 12;

            idx_start = year_offset + m_start;
            idx_end   = year_offset + m_end;

            mask_gs(idx_start : idx_end) = true;
        end

        % Step C
        ndvi1      = v_ndvi(mask_gs);
        vpd_day1   = v_vpdd(mask_gs);
        vpd_night1 = v_vpdn(mask_gs);
        t_day1     = v_td(mask_gs);
        t_night1   = v_tn(mask_gs);
        swc1       = v_swc(mask_gs);
        cld1       = v_cld(mask_gs);

        if length(ndvi1) < 10 || any(isnan(ndvi1))

            valid_idx = ~isnan(ndvi1) & ~isnan(vpd_day1) & ~isnan(vpd_night1) & ...
                ~isnan(t_day1) & ~isnan(t_night1) & ~isnan(swc1) & ~isnan(cld1);

            ndvi1 = ndvi1(valid_idx);
            vpd_day1 = vpd_day1(valid_idx);
            vpd_night1 = vpd_night1(valid_idx);
            t_day1 = t_day1(valid_idx);
            t_night1 = t_night1(valid_idx);
            swc1 = swc1(valid_idx);
            cld1 = cld1(valid_idx);

            if length(ndvi1) < 10; continue; end
        end

        %  Step D
        X = [zscore(detrend(vpd_day1)), ...
            zscore(detrend(vpd_night1)), ...
            zscore(detrend(t_day1)), ...
            zscore(detrend(t_night1)), ...
            zscore(detrend(swc1)), ...
            zscore(detrend(cld1))];
        Y = zscore(detrend(ndvi1));
        try
            [cof, r,pvl,k] = ridge_regress_function(Y, X, k_param, alpha, ridge_k);
            if pvl<0.05 
                row_day_coef(j)   = cof(1);
                row_night_coef(j) = cof(2);
            end
        catch
            row_day_coef(j)   = NaN;
            row_night_coef(j) = NaN;
        end
    end

    ridge_vpd_ndvi_day(i, :)   = row_day_coef;
    ridge_vpd_ndvi_night(i, :) = row_night_coef;

    if mod(i, 10) == 0
        fprintf('Processing row %d / %d\n', i, rows);
    end
end
