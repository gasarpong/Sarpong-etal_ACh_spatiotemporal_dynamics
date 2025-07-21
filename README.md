Sarpong-etal_ACh_spatiotemporal_dynamics

Last update: 2025/7/21

## This README file contains instructions for installing and running the following programs. Example data files can be found in the demo_data folder.
1.	Hierarchical Clustering
2.	3D Spatial cluster mapping and visualization
3.	ROI_temporal_profiles (Peak amplitude, peak latency, dip amplitude, dip latency)
4.	Behavioral Data Analysis: Ymaze_Reward_noRewardTrialResults
5.	Behavioral Data Analysis: findVelocityAndLickRate


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## Hierarchical Clustering
Program: hierarchical_clustering_algorithm.m

This program enables clustering of imaged ROIs using a hierarchical clustering algorithm. 
In brief, clustering is performed following the steps below:

Step 1: Dimensionality reduction is performed using Principal Component Analysis (PCA)

Step 2: Hierarchical clustering is applied on dimensionally reduced data.

## Instructions for using the code:

a. Ensure you have MATLAB installed.

b. Download and extract the matlab file to your preferred directory.

c. Open MATLAB and set the current folder to the directory where the files are located.

d. Add the directory to your MATLAB path using:

addpath(genpath('your_directory_path'));

## Running the Code

a. Open hierarchical_clustering_algorithm.m in MATLAB.

b. Click Run (or press F5).

Follow the prompts in the command window and select your input data.
Note that the input file must contain only the dataset of the desired time window of analysis.

## Output

The script generates output in the form of figures saved in an output file in the same folder as the input data.

## Troubleshooting

If you encounter issues, ensure all dependencies are in the correct path.

Check for missing toolboxes using:

ver

## User-defined parameters:

a. Line 17: retainedVariance = 90; 
Cumulative Variance threshold to decide retained number of Principal components 

b. Line 63: Number of clusters. Either one of these can be used:
i. algorithm-defined number of clusters (using Elbow method)
Ii. user-defined number of clusters

c. Line 216: Computing distance metric using the 'Euclidean' method.
          Using 'complete' method to compute linkageTree 

d. Line 246: Color Threshold for the dendrogram is computed 

e. Line 251: dendrogram(linkageTree,100,'ColorThreshold', threshold);
Number of clusters in the dendrogram can be set by user

The following message will appear once the analysis is completed: "Analysis Completed"
An example dataset ("All_ROIs_clustering") is included. See "demo_data" folder.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## 3D Spatial cluster mapping and visualization
Program: 3D_Spatial_cluster_mapping_Reversal.m

This MATLAB script performs spatial mapping and visualization of region-of-interest (ROI) data across multiple imaging fields (sheets) from an Excel file. Each ROI is mapped onto a standardized 15×15 grid based on its numeric ID, with plotting based on response types and mean activity (e.g., mean z-score or cluster identity).

## Features

- Automatic mapping of ROIs to anatomical space using metadata from each sheet.
- Color-coded open-circle scatter plots based on response type (`Pause`, `Burstpause`, `Pauseburst`, `Others`, or 'NR increase', 'NR decrease', 'others').
- Both 2D and 3D visualizations per sheet and as merged summary plots.
- Flexible to include flipped coordinate systems (`Flip_AP`, `Flip_ML`), .
- Automatically generates a `Results` folder with all plots and a summary Excel file.

## Requirements

- MATLAB R2019b or later
- No additional toolboxes required (standard MATLAB plotting and I/O functions used)

## Input File Format

Provide a single `.xlsx` file with multiple sheets. Each sheet must contain the following columns:
- `AP`, `ML`: Numeric metadata defining the center coordinates.
- `Flip_AP`, `Flip_ML`: Strings indicating axis flipping (`Yes`/`No`).
- Response columns (e.g., `Increase`, `Decrease`, etc.): ROI numbers.
- Mean z-score columns: Must follow format `Mean z score_[response type]`, e.g., `Mean z score_Increase`.
Each ROI number should be numeric, ranging from 1–225 per sheet, assuming a 15×15 grid.

## How to Use

1. Open MATLAB and navigate to the folder containing the script.
2. Run the script directly in the Command Window or Editor.
3. Use the dialog box to select your `.xlsx` input file.
4. The script will process each sheet, generate plots, and save results to:
   ``` /path/to/your/data/Results/ ```
## Outputs include:
   - 2D and 3D `.fig` files per sheet
   - Merged summary 2D and 3D plots
   - Summary Excel file `Summary_data.xlsx` with metadata from all sheets

## Output Example

- `ROI_Map_[sheetname].fig`: 2D spatial map
- `ROI_3D_Map_[sheetname].fig`: 3D map of signal magnitude
- `Merged_ROI_Map_2D.fig`: Combined 2D visualization
- `Merged_ROI_Map_3D.fig`: Combined 3D visualization
- `Summary_data.xlsx`: All ROI coordinates, response types, and values

## Notes

- ROIs are assigned grid positions relative to their numeric ID and sheet index.
- Flipping metadata is used to mirror coordinates across ML or AP axes as needed.
- Empty or malformed ROI entries are skipped with a warning.
- Physical coordinates are centered using the recorded ML and AP metadata values extracted from the input sheets. 
- The script accounts for potential axis inversions resulting from optical properties of the GRIN lens (e.g., mirroring or flipping across ML and AP axes) and image orientation effects introduced by our two-photon imaging system. 
- Each ROI is visualized as an open circle plotted using scatter and scatter3 functions, with edge colors assigned to each functional category.
- When the analysis is completed, a message box will appear with the following message:''Analysis complete. All figures and Excel file are saved in the Results folder.'').



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## ROI_temporal_profiles
Program: ROI_temporal_profiles.m

This program examines fluorescence activity of each ROI, in a pre-determined time window, and extracts the following:
1. Peak amplitude (max activity)
2. Peak latency (time to peak)
3. Dip amplitude (min activity)
4. Dip latency (time to dip)

Note that the input data must contain only the dataset of the time window to be used for the analysis.

Data from different tabs in the same excel file can be analyzed

An example dataset ("ROI_activity") is included. See "demo_data" folder.



## Behavioral Data Analysis
## 1. Program: Ymaze_Reward_noRewardTrialResults.m

% This program analyzes Y-maze behavioral data and generates following plots to show:
      % 1) 'Reward' type or  'No Reward' type  
      % 2) Trial Duration for each trial
      % 3) Rewarded trial percentage for each session, 
      % 4) Number of Rewarded trials / Total trials  for each session
      % 5) Lick occurrences during Reward_Time ± lickDuration (3s etc) period 
      % 6) Lick occurrences during No_Reward_Time  ± lickDuration (3s etc) period 
      % 7) Lick occurrences during first 5s of each trial (user-defined)
      % 8) RewardRate for reward/noReward trials is computed
      % 9) Lick Rate for each 100ms during Reward_Time  ± lickDuration (3s) is found
 
 % A Single file or multiple files can be analyzed simultaneously 

## 2. Program: findVelocityAndLickRate.m

% This program analyzes behavioral data and generates following plots:
      % 1) Velocity at each position  for 'Reward' type or  'No Reward' type 
      % 2) LickRate at each position for 'Reward' type or  'No Reward' type (LickRate for reward period of 2s shows at the end)
      % 3) Each trial's Velocity/ LickRate results are written in output excel files 
      
 % A Single file or multiple files can be analyzed simultaneously 
 % All results are saved to an Output Excel file 

 An example dataset ("Y-maze_behavioral_data") is included. See "demo_data" folder.
