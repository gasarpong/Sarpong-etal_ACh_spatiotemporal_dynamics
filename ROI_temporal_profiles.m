function ROI_temporal_profiles

% This program examines fluorescence activity of each ROI, in a pre-determined time window, and extracts the following:
% 1. Peak amplitude (max activity)
% 2. Peak latency (time to peak)
% 3. Dip amplitude (min activity)
% 4. Dip latency (time to dip)

% The input data should contain only the dataset of the time window to be used for the analysis

% Data from different tabs in the same excel file can be analyzed

    % Open a dialog box to select an Excel file
    [file, path] = uigetfile('*.xlsx', 'Select an Excel file for analysis');
    
    % Check if the user canceled the file selection
    if isequal(file, 0)
        disp('No file selected. Exiting.');
        return;
    end
    
    % Full file path
    filename = fullfile(path, file);
    disp(['File selected: ', filename]);
    
    % Get all worksheet names
    [~, sheets] = xlsfinfo(filename); % Retrieve sheet names
    if isempty(sheets)
        error('No sheets found in the selected Excel file.');
    end
    
    % Initialize output filename
    output_filename = fullfile(path, 'ROI_profiles_results.xlsx');
    
    % Loop through each worksheet
    for s = 1:length(sheets)
        % Read all data (including headers)
        raw_data = readcell(filename, 'Sheet', sheets{s});
        
        % Extract headers (first row)
        headers = raw_data(1, :);
        
        % Extract numeric data (from second row onward)
        numeric_data = cell2mat(raw_data(2:end, :));
        
        % Separate time and ROI activity
        time = numeric_data(:, 1); % First column is time
        roi_activity = numeric_data(:, 2:end); % Remaining columns are ROI activity
        roi_names = headers(2:end); % Extract ROI IDs from headers (skipping 'Time (s)')
        
        % Remove invalid ROI columns (if any)
        valid_indices = ~cellfun(@isempty, roi_names);
        roi_names = roi_names(valid_indices);
        roi_activity = roi_activity(:, valid_indices);
        
        % Validate data consistency
        if isempty(roi_activity) || size(roi_activity, 1) ~= length(time)
            error('Data inconsistency in sheet: %s', sheets{s});
        end

        % Initialize arrays for results
        num_rois = size(roi_activity, 2);
        min_values = zeros(num_rois, 1);
        min_times = zeros(num_rois, 1);
        max_values = zeros(num_rois, 1);
        max_times = zeros(num_rois, 1);

        % Analyze each ROI
        for i = 1:num_rois
            % Find the minimum value and corresponding time
            [min_values(i), min_idx] = min(roi_activity(:, i));
            min_times(i) = time(min_idx);

            % Find the maximum value and corresponding time
            [max_values(i), max_idx] = max(roi_activity(:, i));
            max_times(i) = time(max_idx);
        end

        % Create a table for this sheet
        sheet_results = table(roi_names(:), min_values, min_times, max_values, max_times, ...
                              'VariableNames', {'ROI', 'MinValue', 'MinTime', 'MaxValue', 'MaxTime'});
        
        % Write the table to a corresponding tab in the output Excel file
        writetable(sheet_results, output_filename, 'Sheet', sheets{s});
        disp(['Results saved for sheet: ', sheets{s}]);
    end
    
    disp(['All results saved to: ', output_filename]);
end
