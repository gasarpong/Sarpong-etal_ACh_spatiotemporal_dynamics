function hierarchical_clustering_algorithm()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 % Written by Kavinda Liyanagama, Okinawa Institute of Science and
 % Technology (2025)

 % This program enables clustering of ROI data using the Hierarchical Clustering algorithm

 % Step 1: Dimensionality reduction is performed using Principal Component Analysis (PCA)
 % Step 2: Hierarchical clustering is applied on dimensionally reduced data.

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 clear All;
 clearAllMemoizedCaches; % This MATLAB function clears caches for all MemoizedFunction objects

 retainedVariance = 90; %  Cumulative Variance threshold to decide retained number of Principal components 

 % Prompt the user to select a data file using a file dialog
 [fileName,filePath] = uigetfile('*.*', 'Select a Data File', 'MultiSelect', 'off');  

 % Check if the user selected files
 if isequal(fileName, 0)
    disp('No files selected. Exiting...');
    return;
 end

 [type1,fileName1,ext1] = fileparts(fileName); sheetName = fileName1 ; % Take file name without extension 
 mkdir(strcat(fullfile(filePath,'Results/',sheetName,'HierarchicalClustering'))); % Create a new folder to save results
 outputFileName = strcat(sheetName ,'-HierarChicalClustered_Results.xlsx') ; % Create an output Excel file to save results

 ROIData = readtable(fileName); % Extracted data from the selected sheet 

 % Arrage data for PCA to have each signal in a different row and each data point in a different column 

 noOfROIs = size(ROIData,2);
 data = ROIData{:,1:noOfROIs };
 dataForPCA = transpose(data);  


 % Perform PCA

 [coeff, score, latent, tsquared, explained] = pca(dataForPCA);

 % Plot cumulative variance

 figure1Name= strcat(sheetName,"- Cumulative Variance");
 h1 =figure('Name',figure1Name,'position',[40,150, 700,500]);

 [nComponents ,cumulativeVariancePC] =  findCumulativeVariance(explained,retainedVariance); % Find retained number of Principal components 
 
 dimReducedData = score(:, 1:nComponents); % Dimensionally reduced data

 pc1= score(:,1); % Data of Principal components 1
 pc2 =score(:,2); % Data of Principal components 2
 pc3 =score(:,3); % Data of Principal components 3


 % Cluster data using  Hierarchical Clustering algorithm 

  [linkageTree,noOfClusters,figureList,nameList] = runHierarchicalClustering(dimReducedData,sheetName);
  

 % Choose number of clusters 
  numClusters= noOfClusters ; % Number of clusters ( Decided using Elbow method) 

  %numClusters = 4; % Number of clusters ( Decided by user ) 


 % Cluster data using  Hierarchical Clustering algorithm for decided number of clusters

  clusterLabels = cluster(linkageTree, 'maxclust', numClusters);
  [clusterNumbers,ROICount] = computeCount(clusterLabels); % Compute number of ROIs in each cluster

  % Define custom RGB colors
  customColors = setCustomColors();

 % Visualize clusters (for 2D data)

 figure2Name = strcat(sheetName,"- Clustered Data-2D");
 h2=figure('Name',figure2Name,'position',[120,150, 1000,1000]);

 hold on;
 for i = 1:numClusters
   g(i)= scatter(pc1(clusterLabels == i), pc2(clusterLabels == i), 36, customColors(i, :), 'filled','DisplayName', sprintf('Cluster %d', i));
 end

 setTitleAxis('Hierarchical Clustering on PCA-Reduced Data','Principal Component 1', 'Principal Component 2');

 f2 = plot([0 0], ylim, 'k--', 'LineWidth', 1,'Color', [0.5 0.5 0.5]); % Axis line 1
 f3 = plot(xlim, [0 0], 'k--', 'LineWidth', 1,'Color', [0.5 0.5 0.5]); % Axis line 2

 legend([g(1:numClusters)],'Location', 'best');

 % Display ROI numbers on 2D plot
  labels = strcat(string(1:noOfROIs)); % ROI labels
  %text(pc1,pc2,labels,'VerticalAlignment','bottom','HorizontalAlignment','right')

  hold off;

  %grid on;


 % Visualize clusters (for 3D data)

 figure3Name = strcat(sheetName,"Clustered Data-3D");
 h3= figure('Name',figure3Name,'position',[160,150, 1000,1000]);
 
 plot3DGraph(dimReducedData,numClusters,clusterLabels,customColors);

 % Display ROI numbers on 3D plot
 %text(pc1,pc2,pc3,labels,'VerticalAlignment','bottom','HorizontalAlignment','right')

 hold off;

 grid on;

 % Write results to an output Excel file

 clusteredData = [transpose(1:noOfROIs) clusterLabels];
 clusteredData = array2table ( clusteredData, 'VariableNames' , {'ROI Number','Cluster Number'})  ;  %  ROI number with its cluster number
 writetable(clusteredData,fullfile(filePath,'/','Results/',sheetName,'/','HierarchicalClustering/', outputFileName),'Sheet','Cluster Analysis' ,'Range','A1');

 ROICountData = [clusterNumbers,ROICount]; %  ROI count data for each cluster
 ROICountData = array2table (ROICountData, 'VariableNames' , {'Cluster Number','ROI Count'})  ;  
 writetable(ROICountData,fullfile(filePath,'/','Results/',sheetName,'/','HierarchicalClustering/', outputFileName),'Sheet','Cluster Analysis' ,'Range','D1');
 
 varianceDataPC = [transpose(1:length(explained)) explained]; % How much variance each PC explains 
 clusteredData = array2table (varianceDataPC, 'VariableNames' , {'PC Number','Explained Variance'})  ;  
 writetable(clusteredData,fullfile(filePath,'/','Results/',sheetName,'/','HierarchicalClustering/', outputFileName),'Sheet','Cluster Analysis' ,'Range','H1');

 cumulativeVarianceData = [transpose(1:length(cumulativeVariancePC)) cumulativeVariancePC]; %  Cumulative Variance
 cumulativeVarianceData = array2table (cumulativeVarianceData, 'VariableNames' , {'PC Number','Cumulative Variance'})  ;  
 writetable(cumulativeVarianceData,fullfile(filePath,'/','Results/',sheetName,'/','HierarchicalClustering/', outputFileName),'Sheet','Cluster Analysis' ,'Range','L1');

 variableData = [retainedVariance nComponents numClusters ]; % Variable data
 variableData = array2table (variableData, 'VariableNames' , {'Retained_Variance','No Of PCs Retained', 'NoOfClusters'})  ;  
 writetable(variableData,fullfile(filePath,'/','Results/',sheetName,'/','HierarchicalClustering/', outputFileName),'Sheet','Cluster Analysis' ,'Range','P1');

 outputPCdata  = array2table (dimReducedData, 'VariableNames' , strcat('PC_',string(1:nComponents)))  ;  %  Data for first 3 Principal Components
 writetable(outputPCdata,fullfile(filePath,'/','Results/',sheetName,'/','HierarchicalClustering/', outputFileName),'Sheet','PCData' ,'Range','A1');

 % Save generated figures

 saveas(h1,strcat(fullfile(filePath,'Results/',sheetName,'/','HierarchicalClustering/'),figure1Name,'.fig'))  
 saveas(h2,strcat(fullfile(filePath,'Results/',sheetName,'/','HierarchicalClustering/'),figure2Name,'.fig')) 
 saveas(h3,strcat(fullfile(filePath,'Results/',sheetName,'/','HierarchicalClustering/'),figure3Name,'.fig')) 


 % Save generated figures in sub-functions

 for fileIndex = 1:numel(figureList)
    fig = figureList(fileIndex); 
    saveas(fig,strcat(fullfile(filePath,'Results/',sheetName,'/','HierarchicalClustering/'),nameList(fileIndex),'.fig')) 
 end

 allLabels = transpose(labels); % All ROI Numbers  


 % Initialize an empty cell array to store Cluster Number's data
 clusteredROINumbers = {};

 % Saves raw data of ROIs belong to each cluster to a separate tab the output Excel file

 clusterLabelArray = strings(1, numClusters); % Initialize a string array

 for selectedCluster = 1:numClusters % Loop through each cluster
    selectedROIData = dataForPCA(clusterLabels ==selectedCluster, :);
    selectedROILabels = allLabels(clusterLabels == selectedCluster, :);
    selectedNumberOfROIs = size(selectedROIData,1);

    clusteredROINumbers{selectedCluster} = selectedROILabels;

    outputROIdata  = array2table (transpose(selectedROIData), 'VariableNames' , strcat("ROI ",selectedROILabels(1:selectedNumberOfROIs)))  ;  %  Raw data of selected ROIs
    writetable(outputROIdata,fullfile(filePath,'/','Results/',sheetName,'/','HierarchicalClustering/', outputFileName),'Sheet',strcat("Cluster ",num2str(selectedCluster) ) ,'Range','A1');
  
    clusterLabelArray(selectedCluster) = ['Cluster ' num2str(selectedCluster)]; % Create Cluster labels like 'Cluster 1', 'Cluster 2', ...

 end

 % Write ROI numbers of each cluster to the output Excel file
 
 % Find the maximum number of data points (length) across all signals
 maxLength = max(cellfun(@length, clusteredROINumbers));
 
 % Initialize a matrix (NaN padded) to store all data
 paddedData = NaN(maxLength, length(clusteredROINumbers));

 % Loop through each signal and pad with NaNs
 for i = 1:length(clusteredROINumbers)
    signalLength = length(clusteredROINumbers{i});
    paddedData(1:signalLength, i) = clusteredROINumbers{i};  % Fill in the data
 end

 writematrix(paddedData,fullfile(filePath,'/','Results/',sheetName,'/','HierarchicalClustering/', outputFileName),'Sheet',"Clustered ROIs" ,'Range','A2');
 writematrix(clusterLabelArray, fullfile(filePath,'/','Results/',sheetName,'/','HierarchicalClustering/', outputFileName), 'Sheet',"Clustered ROIs" ,'Range','A1');


% Plot all signals

transposedData = transpose(data);

figureAll = strcat(sheetName,"- All Signals");
figure('Name',figureAll,'position',[160,150, 1400,1000]);


hold on; % Keep all plots in the same figure
for i = 1:size(transposedData, 1)
    selectedColor = clusterLabels(i);

   plot(1:size(transposedData,2), transposedData(i, :),'Color', customColors(selectedColor, :) );
    
end
for i = 1:numClusters 
   plotT1(i)= plot(nan, nan,'Color', customColors(i, :));  
end


% Add labels and title

 setTitleAxis('All Signals Plotted Together','Data Point Number','Z-Score')
 legend(plotT1, strcat("Cluster ",string (1:numClusters) ),'Location', 'best')

 hold off;

 [aveSignalEachCluster, groupAveFigure] = generateGroupFigures(data,transpose(clusterLabels),customColors) ; % Generate figures for each separate group

 % Write data of average signal for each cluster type 

 aveSignalEachCluster  = array2table (transpose(aveSignalEachCluster), 'VariableNames' , strcat("Cluster ",string(1:numClusters)))  ;  %  Average signal of each cluster
 writetable(aveSignalEachCluster,fullfile(filePath,'/','Results/',sheetName,'/','HierarchicalClustering/', outputFileName),'Sheet',"Ave Signal" ,'Range','A1');

 groupAveFigureName = strcat(sheetName,"- Group-Average Signal");
 saveas(groupAveFigure,strcat(fullfile(filePath,'Results/',sheetName,'/','HierarchicalClustering/'),groupAveFigureName,'.fig')) 

end



function [linkageTree,noOfClusters,figureList,nameList] =  runHierarchicalClustering(reducedData,fileName)
% This functions runs Hierarchical clustering on dimensionally reduced data

% Step 1: Compute Distance Matrix

% Compute pairwise distances (e.g., Euclidean)
 distanceMatrix = pdist(reducedData,"euclidean");


% Step 2: Perform Hierarchical Clustering

% Use 'ward' linkage method (can also try ward,'single','average','complete' etc.)
linkageTree = linkage(distanceMatrix,"complete");

% Step 3: Plot the Dendrogram

treeFigureName1= strcat(fileName,"- LinkageTree");
figure1 =figure('Name',treeFigureName1,'position',[120,150, 1200,1200]);

dendrogram(linkageTree,'ColorThreshold', 'default');
title('Hierarchical Clustering Dendrogram (Dimensionality Reduced Data)');
xlabel('Data Point Index');
ylabel('Distance');


% Plot the Dendrogram for all ROIs
treeFigureName2= strcat(fileName,"- LinkageTree -All");
figure2 =figure('Name',treeFigureName2,'position',[120,150, 1200,1200]);

noOfROIs = size(reducedData,1);
%dendrogram(linkageTree,noOfROIs,'ColorThreshold', 2);
dendrogram(linkageTree,100,'ColorThreshold', 4);
title('Hierarchical Clustering Dendrogram (Dimensionality Reduced Data)');
xlabel('Data Point Index');
ylabel('Distance');


% Step 4: Compute Total Within-Cluster Variance
maxClusters = 10; % Maximum number of clusters to evaluate
%withinClusterSum = zeros(maxClusters, 1); % Store within-cluster sum of squares

wcss = zeros(1,maxClusters);  % Array to store WCSS values

% Loop to compute WCSS for each number of clusters
for k = 1:maxClusters
    % Perform hierarchical clustering (using 'complete' linkage)
    Z = linkage(reducedData, 'complete', 'euclidean');
    
    % Assign points to clusters
    T = cluster(Z, 'maxclust', k);
    
    % Calculate the WCSS (within-cluster sum of squares)
    wcss(k) = calculateWCSS(reducedData, T, k);
end


% Step 5: Plot the Elbow Curve

 elbowFigureName = strcat(fileName,"- Elbow Method");
 figure3=figure('Name',elbowFigureName,'position',[80,150, 700,500]);
 hold on;

 noOfClusters= findElbow(1:maxClusters,transpose(wcss) );

 % Pass parameters to save figures 

 figureList = [figure1 figure2 figure3];
 nameList = [treeFigureName1 treeFigureName2 elbowFigureName];

end 



function [numComponents ,cumulativeVariancePC] = findCumulativeVariance(explained,retainedVariance )

% This function computes the explained variance of each principal component, and 
% cumulative variance is computed to decide the number of components needed
% to reach specific cumulative variance (90% etc) decided by the user

% Calculate cumulative variance
 cumulativeVariancePC = cumsum(explained);

% Find the number of components needed for selected cumulative variance ( 90% etc )
 numComponents = find(cumulativeVariancePC >= retainedVariance, 1);

 plot(cumulativeVariancePC(1:15), '-ok','LineWidth', 2);
 hold on;

 setTitleAxis('PCA: Variance Explained','Number of Principal Components', 'Cumulative Variance Explained (%)')
 
 legendString = strcat ('Cumulative : ', num2str(retainedVariance), '%');
 plot(numComponents,cumulativeVariancePC(numComponents), 'ro', 'MarkerSize', 10, 'LineWidth', 2);

 legend('Curve', legendString);
 xticks(1:1:15)
 grid on;
 hold off;

end 



function plot3DGraph(data,numClusters,idx,customColors)

% This function generates a 3D plot with colored data points based on clusters

hold on;

for k = 1:numClusters
    clusterPoints = data(idx == k, :);
    m(k) = scatter3(clusterPoints(:, 1), clusterPoints(:, 2), clusterPoints(:, 3), ...
        50, customColors(k, :), 'o', 'LineWidth', 1.5, 'MarkerFaceColor', 'none', 'DisplayName', sprintf('Cluster %d', k));
end

% Plot the centroids
%  n1 = scatter3(centroids(:, 1), centroids(:, 2), centroids(:, 3), 200, 'k', 'x', ...
%     'LineWidth', 2, 'DisplayName', 'Centroids');

 xlabel('PC 1','fontweight','bold','FontName' , 'Arial' , 'fontsize' , 30);
 ylabel('PC 2','fontweight','bold','FontName' , 'Arial' , 'fontsize' , 30);
 zlabel('PC 3','fontweight','bold','FontName' , 'Arial' , 'fontsize' , 20);
 title('3D Visualization of Hierarchical Clustering');
 %legend('Location', 'best');

 set(gca,'fontweight','bold','fontsize',12); % Set axis ticks' style 
 set(gca, 'LineWidth', 4.5)

 % Get axis limits
 xLimits = xlim;
 yLimits = ylim;
 zLimits = zlim;

% Plot the axis lines
 n2= plot3(xLimits, [0 0], [0 0], 'k--', 'LineWidth', 1,'Color', [0.5 0.5 0.5]); 
 n3= plot3([0 0], yLimits, [0 0], 'k--', 'LineWidth', 1,'Color', [0.5 0.5 0.5]); 
 n4= plot3([0 0], [0 0], zLimits, 'k--', 'LineWidth', 1,'Color', [0.5 0.5 0.5]); 
 legend([m(1:numClusters)],'Location', 'best');

 %grid on;
 view(3); % Set to 3D view

end



function [elbowIndex] = findElbow(x,y)

% This function has been found from internet to find elbow of curve

  plot(x, y, 'o-', 'LineWidth', 2);

% Define the line connecting the first and last points
  p1 = [x(1), y(1)];
  p2 = [x(end), y(end)];

% Calculate distances of all points from the line
  distances = zeros(length(x), 1);

  for i = 1:length(x)
    p = [x(i), y(i)];
    distances(i) = abs(det([p2-p1; p-p1])) / norm(p2-p1); % Point-line distance formula
 end

% Find the index of the maximum distance (elbow point)
 [~, elbowIndex] = max(distances);

% Highlight the elbow point on the plot
 plot(x(elbowIndex), y(elbowIndex), 'ro', 'MarkerSize', 10, 'LineWidth', 2);
 setTitleAxis('Elbow Method for Optimal k','Number of Clusters (k)', 'Sum of Squared Distances');
 xlabel('Number of Clusters (k)');
 legend('Curve', 'Elbow Point');
 grid on;
 hold off;

 % Display elbow point
 %fprintf('Elbow point is at x = %d, y = %d\n', x(elbowIndex), y(elbowIndex));

end



function [clusterNumbers,dataPointCount] = computeCount(clusterIndexes)

% This function computes data point count in each cluster

% Get unique values in cluster Indexes
   clusterNumbers = unique(clusterIndexes);

% Initialize an array to store data point count
  dataPointCount = zeros(size(clusterNumbers));

% Count occurrences of each cluster number
  for i = 1:length(dataPointCount)
    dataPointCount(i) = sum(clusterIndexes == clusterNumbers(i));
  end

end



function wcss = calculateWCSS(data, labels, k)
 % Function to calculate Within-Cluster Sum of Squares(WCSS) for a given clustering solution

    wcss = 0;
    for i = 1:k
        clusterData = data(labels == i, :);  % Get data points in cluster i
        clusterCenter = mean(clusterData, 1);  % Centroid of the cluster
        wcss = wcss + sum(sum((clusterData - clusterCenter).^2));  % Sum of squared distances to centroid
    end
end



function setTitleAxis(plotTitle,xTitle, yTitle)
   % This function customizes plot view
   
    title(plotTitle)   
    ylabel(yTitle,'fontweight','bold','FontName' , 'Arial' , 'fontsize' , 12)
    xlabel(xTitle,'fontweight','bold','FontName' , 'Arial' , 'fontsize' , 12)  
      
    ax = gca;
    set(gca,'fontweight','bold','fontsize',12); % Set axis ticks' style  
end



function [aveSignalEachGroup, groupAveFigure]= generateGroupFigures(data,groupType,customColors)

% This function generates figures with plots for each cluster type separately

% Parameters
maxSubplots = 30; % Maximum subplots per figure (5x4 grid)
figureWidth = 1300; % Width of the figure
figureHeight = 1000; % Height of the figure
uniqueGroups = unique(groupType); % Unique group types
maxPlots = 20 ;  % Set a maximum plot limit for each cluster type to prevent program crashes caused by generating too many figures.

aveSignalEachGroup =[];
ROICountEachGroup =[];

% Loop through each group
for group = uniqueGroups

    groupCols = find(groupType == group); % Find columns belonging to this group
    numColsInGroup = length(groupCols); % Number of columns in this group

    selectedDataAll = data(:,groupCols);
    aveSignal = mean(selectedDataAll, 2);
    aveSignalEachGroup(group,:) = aveSignal;
    ROICountEachGroup(group) = numColsInGroup;

    if (numColsInGroup > maxPlots)
        numColsInGroup =maxPlots; % Set a maximum plot limit for each cluster type to prevent program crashes caused by generating too many figures.
    end
    
    % Iterate over columns in the group
    for idx = 1:numColsInGroup
        col = groupCols(idx); % Current column index in data

        % Determine figure and subplot position
        figureIndex = ceil(idx / maxSubplots); % Current figure number for this group
        subplotIndex = mod(idx-1, maxSubplots) + 1; % Subplot position in the figure
        
        % Open a new figure when starting a new batch of subplots
        if subplotIndex == 1
            figure % Create new figure
            %figureName = 'Group '+ num2str(group) +' - Figure '+ num2str(figureIndex);

            set(gcf, 'Position', [100, 100, figureWidth, figureHeight]); % Set figure size
            sgtitle(['Cluster ' num2str(group) ' - Figure ' num2str(figureIndex)]); % Super title
            pause(1);  % Pauses execution for 1 second to prevent crashing due to too many figures 
        end
        
        % Add subplot to the current figure
          subplot(6, 5, subplotIndex); % 5x4 grid
          plot(data(:, col),'Color', customColors(group, :)); % Plot the column data
          xlim([0 size(data(:, col),1)]);
        
          title("ROI "+ num2str(col)); % Title for each subplot
           % xlabel('X-axis');
           % ylabel('Y-axis');
    end 
end

% Plot grouped-average signals
figureGroupAve = strcat("Group-Average Signals");
groupAveFigure =figure('Name',figureGroupAve,'position',[160,150, 800,1000]);

hold on; % Keep all plots in the same figure

noOfClusters = size(aveSignalEachGroup, 1);
    
for i = 1:noOfClusters 
   plotT1(i)= plot(1:size(aveSignalEachGroup,2), aveSignalEachGroup(i, :),'Color', customColors(i, :)); 
    
end

for i = 1:numel(plotT1)
    %Add a new row to DataTip to show the DisplayName of the line
     % plotT1(i).DataTipTemplate.DataTipRows(i) = dataTipTextRow( 'Cluster ',num2str(i));  

    plotT1(i).DataTipTemplate.DataTipRows(i) = dataTipTextRow('Cluster' ,repmat({string(i)},size( plotT1(i).XData)));    
end 

setTitleAxis('Average Signal for each group','Data point number', 'Z-Score');


% Customize legend names
legend( strcat("Cluster ", string (1:size(aveSignalEachGroup, 1)), "  - (",string(ROICountEachGroup(1:size(aveSignalEachGroup, 1))), ")" ), 'Location', 'best');    
hold off;

msgbox('Analysis Completed!', 'Alert');

end


function [customColors] = setCustomColors()
% This function defines custom RGB colors

customColors = [
    1, 0, 1;        % Magenta
    0, 1, 0;        % Green
    0, 0, 1;        % Blue
    0, 0, 0;        % Black
    1, 0.65, 0;     % Orange
    0, 1, 1;        % Cyan
    0, 0, 1;        % Blue
    1, 0, 0;        % Red
    1, 0.65, 0;     % Orange
    0.5, 0.5, 0.5;  % Gray
    0.5, 0, 0;      % Dark Red
    0, 0.5, 0;      % Dark Green
    0.75, 0, 0.75;  % Purple
    0.5, 0.5, 1;    % Light Blue
    0.25, 0.25, 0;  % Olive
    0.5, 1, 0.5;    % Light Green
    1, 1, 0;        % Yellow
    0.75, 0.75, 0;  % Light Yellow
    0.2, 0.6, 1;    % Sky Blue
    0.8, 0.8, 0.8;  % Light Gray
    0.4, 0.2, 0.6;  % Lavender
    0.7, 0.3, 0.2;  % Brown
 ];

end
