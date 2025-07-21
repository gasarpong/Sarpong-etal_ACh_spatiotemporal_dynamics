% Load Excel file with multiple sheets
[filename, pathname] = uigetfile({'*.xlsx;*.xls'}, 'Select the Excel File');
if isequal(filename, 0), disp('User canceled.'); return; end
filePath = fullfile(pathname, filename);
sheets = sheetnames(filePath);

% Create Results folder
resultsDir = fullfile(pathname, 'Results');
if ~exist(resultsDir, 'dir')
    mkdir(resultsDir);
end

% ROI grid parameters
grid_size = [15, 15];
ROI_spacing = 50 * 0.001;
shrink_factor = 1;

% Define colors
colors = struct('Increase', [1 0 0], 'Decrease', [0 0 1], 'Nonresponsive', [1 0.5 0]);
edge_colors = struct('Increase', [1.0 0.5 0.0], 'Decrease', [0 0 1], 'Nonresponsive', [0.5 1 0.5]);

marker_size = 50;
marker_linewidth = 1;

merged_data = [];

h = waitbar(0, 'Processing data...');
total_sheets = length(sheets);

for s = 1:total_sheets
    waitbar(s/total_sheets, h, ['Processing sheet: ' sheets{s}]);

    data = readtable(filePath, 'Sheet', sheets{s}, 'PreserveVariableNames', true);
    AP_center = unique(data.AP); ML_center = unique(data.ML);
    Flip_AP = string(unique(data.Flip_AP)); Flip_ML = string(unique(data.Flip_ML));

    [X, Y] = meshgrid(linspace(-7, 7, grid_size(1)), linspace(-7, 7, grid_size(2)));
    X = X * ROI_spacing * shrink_factor + ML_center;
    Y = Y * ROI_spacing * shrink_factor + AP_center;

    if strcmpi(Flip_ML, 'Yes'), X = 2 * ML_center - X; end
    if strcmpi(Flip_AP, 'Yes'), Y = 2 * AP_center - Y; end

    sheet_data = [];
    z_values = []; x_values = []; y_values = [];
    c_values = []; edge_list = {};
    total_count = 0;
    type_counts = struct('Increase', 0, 'Decrease', 0, 'Nonresponsive', 0);

    figure('Visible', 'on'); hold on;
    title(['ROI Map: ', sheets{s}]);
    xlabel('Mediolateral (mm)'); ylabel('Anteroposterior (mm)');
    set(gca, 'FontSize', 12, 'FontName', 'Arial', 'LineWidth', 1.5);

    for type = {'Increase', 'Decrease', 'Nonresponsive'}
        type_str = type{1};
        mean_z_col = ['Mean z score_' type_str];

        if ~ismember(type_str, data.Properties.VariableNames) || ...
           ~ismember(mean_z_col, data.Properties.VariableNames)
            continue;
        end

        ROI_list = data.(type_str);
        z_scores = data.(mean_z_col);

        if ~iscell(ROI_list), ROI_list = cellstr(string(ROI_list)); end
        if ~iscell(z_scores), z_scores = num2cell(z_scores); end

        valid_idx = ~cellfun(@isempty, ROI_list) & ~cellfun(@isempty, z_scores);
        ROI_list = ROI_list(valid_idx);
        z_scores = cell2mat(z_scores(valid_idx));

        type_counts.(type_str) = length(ROI_list);
        total_count = total_count + length(ROI_list);

        for j = 1:length(ROI_list)
            roi_str = lower(strtrim(ROI_list{j}));
            roi_num = regexp(roi_str, '\d+', 'match');

            if ~isempty(roi_num)
                roi_num = str2double(roi_num{1});
                if roi_num <= prod(grid_size)
                    row = floor((roi_num - 1) / grid_size(1)) + 1;
                    col = mod((roi_num - 1), grid_size(2)) + 1;

                    this_color = colors.(type_str);
                    this_edge = edge_colors.(type_str);
                    this_z = z_scores(j);

                    z_values = [z_values; this_z];
                    x_values = [x_values; X(row, col)];
                    y_values = [y_values; Y(row, col)];
                    c_values = [c_values; this_color];
                    edge_list{end+1} = this_edge;

                    scatter(X(row, col), Y(row, col), marker_size, this_color, ...
                        'o', 'MarkerEdgeColor', this_edge, 'LineWidth', marker_linewidth);

                    sheet_data = [sheet_data; {sheets{s}, roi_num, X(row, col), ...
                        Y(row, col), type_str, this_z, this_color}];
                end
            end
        end
    end

    legend_labels = {};
    for type = {'Increase', 'Decrease', 'Nonresponsive'}
        type_str = type{1};
        pct = 100 * type_counts.(type_str) / max(total_count, 1);
        legend_labels{end+1} = sprintf('%s (%.1f%%)', type_str, pct);
    end
    xlim([floor(min(x_values)*10)/10, ceil(max(x_values)*10)/10]);
    ylim([floor(min(y_values)*10)/10, ceil(max(y_values)*10)/10]);
    legend(legend_labels, 'Location', 'northeastoutside');
    hold off;
    saveas(gcf, fullfile(resultsDir, ['ROI_Map_' sheets{s} '.fig']));

    % 3D Plot
    figure('Visible', 'on'); hold on;
    %title(['3D ROI Activity Map: ', sheets{s}]);
    %xlabel('Mediolateral (mm)'); ylabel('Anteroposterior (mm)'); zlabel('Mean ACh (z)');
    set(gca, 'FontSize', 30, 'FontName', 'Arial', 'LineWidth', 3.5);
    for i = 1:length(z_values)
        scatter3(x_values(i), y_values(i), z_values(i), marker_size, ...
            c_values(i,:), 'o', 'MarkerEdgeColor', edge_list{i}, 'LineWidth', marker_linewidth);
    end
    view(45, 30); grid on;
    saveas(gcf, fullfile(resultsDir, ['ROI_3D_Map_' sheets{s} '.fig']));
    hold off;

    merged_data = [merged_data; sheet_data];
end

% ===== Merged Plotting and Excel Export =====
% Convert to cell if table
if istable(merged_data)
    sorted_data = table2cell(merged_data);
else
    sorted_data = merged_data;
end

% Sort by type: Increase > Decrease > Nonresponsive
priority = containers.Map({'Increase', 'Decrease', 'Nonresponsive'}, [1, 2, 3]);
type_order = cellfun(@(t) priority(t), sorted_data(:,5));
[~, sort_idx] = sort(type_order);
sorted_data = sorted_data(sort_idx,:);

% 2D Merged Plot
figure('Visible', 'on'); hold on;
title('Merged ROI Map');
xlabel('Mediolateral (mm)'); ylabel('Anteroposterior (mm)');
set(gca, 'FontSize', 12, 'LineWidth', 1.5);
for i = 1:size(sorted_data,1)
    type_str = sorted_data{i,5};
    edge_color = edge_colors.(type_str);
    scatter(sorted_data{i,3}, sorted_data{i,4}, marker_size, ...
        'MarkerEdgeColor', edge_color, 'LineWidth', marker_linewidth);
end
legend(legend_labels, 'Location', 'northeastoutside');
saveas(gcf, fullfile(resultsDir, 'Merged_ROI_Map_2D.fig'));
hold off;

% 3D Merged Plot
figure('Visible', 'on'); hold on;
%title('Merged 3D ROI Activity Map');
%xlabel('Mediolateral (mm)'); ylabel('Anteroposterior (mm)'); zlabel('Mean ACh (z)');
set(gca, 'FontSize', 30, 'LineWidth', 3.5);
for i = 1:size(sorted_data,1)
    x = sorted_data{i,3};
    y = sorted_data{i,4};
    z = sorted_data{i,6};
    type_str = sorted_data{i,5};
    edge_color = edge_colors.(type_str);
    scatter3(x, y, z, marker_size, 'o', 'MarkerEdgeColor', edge_color, ...
        'LineWidth', marker_linewidth);
end
view(45, 30); grid on;
saveas(gcf, fullfile(resultsDir, 'Merged_ROI_Map_3D.fig'));
hold off;

% ===== Save to Excel: 
% Reorder columns: 'Sheet', 'ROI', 'Type', 'X', 'Y', 'Z', 'ColorRGB'
reordered_data = cell(size(sorted_data,1), 7);
for i = 1:size(sorted_data,1)
    reordered_data(i,:) = { ...
        sorted_data{i,1}, ... % Sheet
        sorted_data{i,2}, ... % ROI
        sorted_data{i,5}, ... % Type
        sorted_data{i,3}, ... % X
        sorted_data{i,4}, ... % Y
        sorted_data{i,6}, ... % Z
        sorted_data{i,7}      % ColorRGB
    };
end

headers = {'Sheet', 'ROI', 'Type', 'X', 'Y', 'Z', 'ColorRGB'};
output_excel = fullfile(resultsDir, 'Summary_data.xlsx');

% Write merged data to the first sheet
writecell([headers; reordered_data], output_excel, 'Sheet', 'All_Merged');

% Write individual sheet data to corresponding Excel sheets
unique_sheets = unique(sorted_data(:,1));
for i = 1:length(unique_sheets)
    sheet_name = unique_sheets{i};
    this_data = reordered_data(strcmp(sorted_data(:,1), sheet_name), :);
    writecell([headers; this_data], output_excel, 'Sheet', sheet_name);
end

% Final message
close(h);
msgbox('Analysis complete. All figures and Excel file are saved in the Results folder.', 'Done');