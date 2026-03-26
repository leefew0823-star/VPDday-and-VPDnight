clear; clc;
input_file = '---';
var_name   = '---';
start_year = 1982;
end_year   = 2022;
n_years    = end_year - start_year + 1; 
obs_per_yr = 24;

%% SG
sg_order = 3;   
sg_window = 5;  

%% Dynamic Thresholds
sos_thresh = 0.2; % 20% amplitude for Start
eos_thresh = 0.5; % 50% amplitude for End
amp_min    = 0.1; 

%% Load Data

data_struct = load(input_file);
if isfield(data_struct, var_name)
    ndvi_data = data_struct.(var_name);
else
    vars = fieldnames(data_struct);
    ndvi_data = data_struct.(vars{1});
end

[rows, cols, time_steps] = size(ndvi_data);

SOS_Matrix = NaN(rows, cols, n_years);
EOS_Matrix = NaN(rows, cols, n_years);

%% Data Processing

doy_15 = round(linspace(1, 365, obs_per_yr)); 
doy_daily = 1:365;

poolobj = gcp('nocreate');
if isempty(poolobj)
    parpool; 
end

parfor r = 1:rows
    row_sos = NaN(cols, n_years);
    row_eos = NaN(cols, n_years);
    
    for c = 1:cols
        % Step A: 
        ts_raw = double(squeeze(ndvi_data(r, c, :)));

        if all(isnan(ts_raw)) || mean(ts_raw, 'omitnan') < 0.05
            continue; 
        end
        
        %  Step B:
        ts_filled = fillmissing(ts_raw, 'linear');
        
        % Step C
        try
            ts_smooth = sgolayfilt(ts_filled, sg_order, sg_window);
        catch
            continue;
        end
        
        % Step D:
        for y = 1:n_years
            
            idx_start = (y-1) * obs_per_yr + 1;
            idx_end   = idx_start + obs_per_yr - 1;
            
            
            yr_data_15 = ts_smooth(idx_start:idx_end);
            
            % Step E
            yr_data_daily = interp1(doy_15, yr_data_15, doy_daily, 'spline');
            
            % Step F
            
            [pks, loc_max] = max(yr_data_daily);
            
           
            min_val = min(yr_data_daily);
            amp = pks - min_val;

            if amp < amp_min
                continue;
            end
            
           
            thresh_sos_val = min_val + sos_thresh * amp;
            
            left_curve = yr_data_daily(1:loc_max);
            
            idx_sos = find(left_curve >= thresh_sos_val, 1, 'first');
            
            if ~isempty(idx_sos)
                row_sos(c, y) = idx_sos;
            end
            
            
            thresh_eos_val = min_val + eos_thresh * amp;
            
            right_curve = yr_data_daily(loc_max:end);
            
            idx_eos_local = find(right_curve >= thresh_eos_val, 1, 'last');
            
            if ~isempty(idx_eos_local)
                row_eos(c, y) = loc_max + idx_eos_local - 1;
            end
        end
    end
    
   
    SOS_Matrix(r, :, :) = row_sos;
    EOS_Matrix(r, :, :) = row_eos;

    if mod(r, 10) == 0
        fprintf('Processing Row %d / %d...\n', r, rows);
    end
end

%% Save Results
outfile_sos = sprintf('---.mat', start_year, end_year);
outfile_eos = sprintf('---.mat', start_year, end_year);
save(outfile_sos, 'SOS_Matrix', '-v7.3'); 
save(outfile_eos, 'EOS_Matrix', '-v7.3');
