 function Ymaze_findVelocityAndLickRate()
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
 
 % Written by Kavinda Liyanagama, Okinawa Institute of Science and Technology (2025)

 % This program analyzes data and generates the following plots 
 
      % 1) Velocity at each position  for 'Reward' type or  'No Reward' type , 
      % 2) LickRate at each position for 'Reward' type or  'No Reward' type (LickRate for reward period of 2s shows at the end )
      % 3) Each trial's Velocity/ LickRate results are written output excel files 

      
 % A Single file or multiple simultaneous files can be analyzed with ths program.
 % All results are saved to output Excel files.
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  close All
  clearAllMemoizedCaches; % This MATLAB function clears caches for all MemoizedFunction objects
  
  noOfPositions = 9; % No of positions in each Trial
  
  %  Open file and take files' name and path  
  [files,filePath] = uigetfile('*.*','MultiSelect','on'); % Select multiple files 
  if ~iscell(files)
      files = {files}; % Convet to a cell array if only 1 file is selected 
  end
 
  noOfSelectedFiles= length(files); % Find number of selected files
    
  if  noOfSelectedFiles ==0 % Check if no Behavioural file is selected
     return;      % Quit the function 
  end
  
  
  [typeFile1,fileName1,~] = fileparts(files{1}); % Take the first file name without extension 
  [typeFile2,fileName2,~] = fileparts(files{noOfSelectedFiles}); % Take the last file name without extension 
  joinedFileName = strcat(fileName1,'-',fileName2);
    
  if noOfSelectedFiles == 1                
       joinedFileName = fileName1;   % If only one file is selected 
  end 
   
  mkdir(strcat(fullfile(filePath,'Results/',joinedFileName))); % Create a new folder to save results
  outputFileName1 = strcat(joinedFileName ,'- Velocity_LickRateResults.xls') ;  % Create an output Excel file to save results  
  outputFileName2 = strcat(joinedFileName ,'- Velocity_Reward_All.xls') ; 
  outputFileName3 = strcat(joinedFileName ,'- Velocity_NoReward_All.xls') ; 
  outputFileName4 = strcat(joinedFileName ,'- LickRate_Reward_All.xls') ; 
  outputFileName5 = strcat(joinedFileName ,'- LickRate_NoReward_All.xls') ; 
  
    
  x0=100; y0=100; width=600;height=400; % Set origin and dimensions of figures
      
  for i=1:noOfSelectedFiles  %  Select each file 

      fileName = files{i};
      [type1,fileNameShort,~] = fileparts(fileName); % Take file name without extension 
      
      [~,SheetNames] = xlsfinfo(fileName); % Find sheet names
         
      opts = detectImportOptions(fileName); 
      opts.SelectedVariableNames= {'DateTime','sense1Events','SystemMsg','MsgValue1','MsgValue2'}; % Read data from these selected columns
      filteredData = readtable(fileName, opts); 
      [noOfRows,noOfColums] = size(filteredData);   % Filtered  array size 

      eventData = filteredData.(3);  % Extract SystemMsg data 
      eventNameData = eventData( strcmp(eventData,'Reward')| strcmp(eventData,'No Reward') ==1,:); % Extract rows that have 'Reward' or ' No Reward' string 
 
      Index_Reward   = transpose(find(~contains(eventNameData,'No Reward'))); % Find indexes of 'Reward' trials 
      Index_noReward = transpose(find(contains(eventNameData,'No Reward')));  % Find indexes of 'No Reward' trials
      
      velocityRewardArray = zeros(length(Index_Reward),noOfPositions-1);  % Array to store Velocity at each position for all Reward trials 
      velocityNoRewardArray = zeros(length(Index_noReward), noOfPositions-1); % Array to store Velocity at each position for all NoReward trials
      
      lickRateRewardArray = zeros(length(Index_Reward),noOfPositions); % Array to store LickRate at each position for all Reward trials 
      lickRateNoRewardArray = zeros(length(Index_noReward), noOfPositions); % Array to store LickRate at each position for all NoReward trials
      
      currentRowNo= 1;  % Start to read data from this row
      rewardTrialNumber = 1;  % Initial Reward trial number
      noRewardTrialNumber =1; % Initial No Reward trial number
      
      positionArray = (40:40:320) ; % Position Array ( startPosition: positionDistance : endPosition )
      positionArray2 =(40:40:360) ;  % Position Array2 (last position is allocated for post- reward position licks   )
     
      while currentRowNo <= noOfRows % Loop through each row
              
        if contains(filteredData.(3)(currentRowNo),'start trial')==1  % Check if it's a trial start
            
            currentRowNo = currentRowNo+1; % Next row ( Position 0)
                                                                        % findLickingRateVelocity(currentRowNo,filteredData,noOfPositions);
           [tempVelocityArray,tempLickRateArray,trialType,currentRowNo]=  findLickingRateVelocity(currentRowNo,filteredData,noOfPositions);
           
           if trialType == 1 % Check of it's a NoReward Trial
               velocityNoRewardArray(noRewardTrialNumber,:)= tempVelocityArray;
               lickRateNoRewardArray(noRewardTrialNumber,:)= tempLickRateArray;
               noRewardTrialNumber = noRewardTrialNumber + 1;
               
           elseif trialType == 2 % Check of it's a Reward Trial
               velocityRewardArray(rewardTrialNumber,:)= tempVelocityArray;
               lickRateRewardArray(rewardTrialNumber,:)= tempLickRateArray;
               rewardTrialNumber = rewardTrialNumber + 1;
               
           end
           
        end   
           
        currentRowNo=currentRowNo+1; % Next row
    
      end
      
      [meanVelocityRewardArray,velocityReward_SEM] = findMeanSEM(velocityRewardArray);
      [meanVelocityNoRewardArray,velocityNoReward_SEM] = findMeanSEM(velocityNoRewardArray);
      
      [meanLickRateRewardArray,lickRateReward_SEM] = findMeanSEM(lickRateRewardArray);
      [meanLickRateNoRewardArray,lickRateNoReward_SEM] = findMeanSEM(lickRateNoRewardArray); 
      
      
      figureNameVelocity = [fileNameShort,'- Velocity Results']; 
      x0 = x0+20;  y0= y0+20;
      h1 = figure('Name',figureNameVelocity,'position',[x0,y0,width,height]);
      
      errorbar( meanVelocityRewardArray,positionArray,velocityReward_SEM,'horizontal' , 'g-o', 'MarkerSize',8 ,'MarkerFaceColor','green' ,'MarkerEdgeColor','green'); hold on 
      errorbar( meanVelocityNoRewardArray,positionArray,velocityNoReward_SEM,'horizontal', 'k-o','MarkerSize',8, 'MarkerFaceColor','black','MarkerEdgeColor','black'); hold on
   
      legend('Reward','NoReward','Location', 'Best'); 
     
      title('Reward/NoReward - Velocity')   
      xlabel('Velocity ( cm/s )')  
      ylabel('Position ( cm )')  
      set(gca,'YTick',0:40:320,'fontweight','bold','fontsize',12)
      % set(gca,'LooseInset',get(gca,'TightInset')); % Set figure's outerspace 
      ylim([1 350]);
      
      figureNameLickRate = [fileNameShort,'- LickRate Results']; 
      h2 = figure('Name',figureNameLickRate,'position',[x0+700,y0,width,height]);
      
      errorbar(meanLickRateRewardArray,positionArray2,lickRateReward_SEM,'horizontal' ,'r-o', 'MarkerSize',8 ,'MarkerFaceColor','red' ,'MarkerEdgeColor','red'); hold on 
      errorbar(meanLickRateNoRewardArray,positionArray2,lickRateNoReward_SEM,'horizontal' , 'k-o','MarkerSize',8, 'MarkerFaceColor','black','MarkerEdgeColor','black'); hold off
   
      legend('Reward','NoReward','Location', 'Best')
      title('Reward/NoReward - LickRate')   
      xlabel('LickRate ( No Of Licks /s )')  
      ylabel('Position ( cm )')  
      yline( 320, '--b','HandleVisibility','off' ); % Mark  reward position 
      y_labels = [40:40:320, "Outcome"]; % Add "Outcome" label to Reward period
      
      set(gca,'YTick',40:40:360,'YTickLabel',y_labels,'fontweight','bold','fontsize',12)
      % set(gca,'LooseInset',get(gca,'TightInset')); % Set figure's outerspace 
      ylim([1 370]);
      
      outputDataArray1 = [transpose(positionArray) transpose(meanVelocityRewardArray) transpose(meanVelocityNoRewardArray) ] ; % Take all results 
      outputDataTable1 = array2table (outputDataArray1, 'VariableNames' , {'Positions' , 'Mean_Velocity_Reward', 'Mean_Velocity_NoReward'})  ;  % Convert array to table    
      writetable(outputDataTable1,fullfile(filePath,'Results/',joinedFileName,'/',outputFileName1),'sheet',SheetNames{1}  ,'Range','A1') % Write results to an output excel file   
      
      outputDataArray2 = [transpose(positionArray2)  transpose(meanLickRateRewardArray) transpose(meanLickRateNoRewardArray) ] ; % Take all results 
      outputDataTable2 = array2table (outputDataArray2, 'VariableNames' , {'Positions' , 'Mean_LickRate_Reward', 'Mean_LickRate_NoReward'})  ;  % Convert array to table    
      writetable(outputDataTable2,fullfile(filePath,'Results/',joinedFileName,'/',outputFileName1),'sheet',SheetNames{1}  ,'Range','E1') % Write results to the output excel file   

      outputDataArray3 = [transpose(velocityReward_SEM) transpose(velocityNoReward_SEM) ] ; % Take all results 
      outputDataTable3 = array2table (outputDataArray3, 'VariableNames' , {'SEM_Velocity_Reward', 'SEM_Velocity_NoReward'})  ;  % Convert array to table    
      writetable(outputDataTable3,fullfile(filePath,'Results/',joinedFileName,'/',outputFileName1),'sheet',SheetNames{1}  ,'Range','I1') % Write results to the output excel file   

      outputDataArray4 = [transpose(lickRateReward_SEM) transpose(lickRateNoReward_SEM) ] ; % Take all results 
      outputDataTable4 = array2table (outputDataArray4, 'VariableNames' , {'SEM_LickRate_Reward','SEM_LickRate_NoReward'})  ;  % Convert array to table    
      writetable(outputDataTable4,fullfile(filePath,'Results/',joinedFileName,'/',outputFileName1),'sheet',SheetNames{1}  ,'Range','L1') % Write results to an output excel file   
      writematrix("Outcome",fullfile(filePath,'Results/',joinedFileName,'/',outputFileName1),'sheet',SheetNames{1}  ,'Range','E10') % Write "Outcome" to the output excel file   
 
      writematrix("Position",fullfile(filePath,'Results/',joinedFileName,'/',outputFileName2),'sheet',SheetNames{1}  ,'Range','A1') 
      writematrix(transpose(positionArray),fullfile(filePath,'Results/',joinedFileName,'/',outputFileName2),'sheet', SheetNames{1} ,'Range','A2');   
      outputDataTable5 = array2table (transpose(velocityRewardArray), 'VariableNames' , strcat('Trial_',string(Index_Reward(1:end))))  ;  %  Velocity for each Reward trial
      writetable(outputDataTable5,fullfile(filePath,'Results/',joinedFileName,'/',outputFileName2),'sheet', SheetNames{1} ,'Range','B1');  
      
      writematrix("Position",fullfile(filePath,'Results/',joinedFileName,'/',outputFileName3),'sheet',SheetNames{1}  ,'Range','A1') 
      writematrix(transpose(positionArray),fullfile(filePath,'Results/',joinedFileName,'/',outputFileName3),'sheet', SheetNames{1} ,'Range','A2');   
      outputDataTable6 = array2table (transpose(velocityNoRewardArray), 'VariableNames' , strcat('Trial_',string(Index_noReward(1:end))))  ;  %  Velocity for each No-Reward trial
      writetable(outputDataTable6,fullfile(filePath,'Results/',joinedFileName,'/',outputFileName3),'sheet', SheetNames{1} ,'Range','B1');                      

      writematrix("Position",fullfile(filePath,'Results/',joinedFileName,'/',outputFileName4),'sheet',SheetNames{1}  ,'Range','A1') 
      writematrix(transpose(positionArray2),fullfile(filePath,'Results/',joinedFileName,'/',outputFileName4),'sheet', SheetNames{1} ,'Range','A2');   
      outputDataTable7 = array2table (transpose(lickRateRewardArray), 'VariableNames' , strcat('Trial_',string(Index_Reward(1:end))))  ;  %  Velocity for each Reward trial
      writetable(outputDataTable7,fullfile(filePath,'Results/',joinedFileName,'/',outputFileName4),'sheet', SheetNames{1} ,'Range','B1');  
      writematrix("Outcome",fullfile(filePath,'Results/',joinedFileName,'/',outputFileName4),'sheet',SheetNames{1}  ,'Range','A10') % Write "Outcome" to the output excel file   
 
      writematrix("Position",fullfile(filePath,'Results/',joinedFileName,'/',outputFileName5),'sheet',SheetNames{1}  ,'Range','A1') 
      writematrix(transpose(positionArray2),fullfile(filePath,'Results/',joinedFileName,'/',outputFileName5),'sheet', SheetNames{1} ,'Range','A2');   
      outputDataTable8 = array2table (transpose(lickRateNoRewardArray), 'VariableNames' , strcat('Trial_',string(Index_noReward(1:end))))  ;  %  Velocity for each No-Reward trial
      writetable(outputDataTable8,fullfile(filePath,'Results/',joinedFileName,'/',outputFileName5),'sheet', SheetNames{1} ,'Range','B1');                      
      writematrix("Outcome",fullfile(filePath,'Results/',joinedFileName,'/',outputFileName5),'sheet',SheetNames{1}  ,'Range','A10') % Write "Outcome" to the output excel file   
 
      saveas(h1,strcat(fullfile(filePath,'Results/',joinedFileName,'/'),figureNameVelocity,'.fig')) % Save generated figure
      saveas(h2,strcat(fullfile(filePath,'Results/',joinedFileName,'/'),figureNameLickRate,'.fig')) % Save generated figure
              
   end 
 end 
 
 
 function [ tempVelocityArray,tempLickRateArray,trialType,nextRowNo]=  findLickingRateVelocity(startRowNo,data,noOfPositions)
    % This function returns 'Velocity info and Licking Rate' info at each position of current Trial
 
    timeStart = data.(1)(startRowNo); % Find current Trial's time at position 0
    tempVelocityArray =  zeros(1,noOfPositions-1); % To store Velocity at each position
    tempLickRateArray =  zeros(1,noOfPositions-1); % To store LickRate at each position
    
    currentRowNo = startRowNo; % To set each Position's row number
    distancePreviousPosition = 0; % Set initial Position value 
    timePreviousPosition = timeStart; % initialtime value
    tempPositionIndex =1; % Index of position 0
    trialType = 0; % Default trial type
    
    while ( ~ (contains(data.(3)(currentRowNo),'start trial')|| contains(data.(3)(currentRowNo),'end')) )  % Continue until next trial or  'end' is found 
        
        if ( contains(data.(4)(currentRowNo),  strcat( 'Position',num2str(tempPositionIndex )))) % Check if current row has 'Position data'
          
           if tempPositionIndex >1  % Check if it's not Position 0
                 
              timeCurrentPosition = data.(1)(currentRowNo); % Time of current position 
              timeDiffPositions = etime(datevec(timeCurrentPosition),datevec(timePreviousPosition)); % Find time to current Position from Previous Position 
                          
              % Find velocity at current position 
              
              distanceDiffPositions =  data.(5)(currentRowNo) - distancePreviousPosition ; % Find distance between current position and previous position
              velocityTemp = distanceDiffPositions / timeDiffPositions; % Find Velocity of mouse between current position and previous position
              tempVelocityArray(tempPositionIndex-1) = velocityTemp;  
             
              % Find LickRate at current position 
           
              noOfLicks = nnz( data.(2) ((etime( datevec(data.(1)), datevec(timePreviousPosition )) <=etime(datevec(timeCurrentPosition) , datevec(timePreviousPosition )))  & (etime( datevec(data.(1)), datevec(timePreviousPosition )) >=0) & (data.(2) ==1)  ));
              lickRateTemp = noOfLicks/timeDiffPositions; % Find lick rate from last position to  current position 
              tempLickRateArray(tempPositionIndex-1) = lickRateTemp;
                  
              distancePreviousPosition =  data.(5)(currentRowNo);
              timePreviousPosition = timeCurrentPosition; 
              
              if tempPositionIndex == noOfPositions % Check if mouse has reached rewarding position
                  
                  % Find lick rate for the period of 2s since mouse reaches reward position 
                  
                  noOfLicksPostReward = nnz( data.(2) ((etime( datevec(data.(1)), datevec(timePreviousPosition )) <= 2 )  & (etime( datevec(data.(1)), datevec(timePreviousPosition )) >=0) & (data.(2) ==1)  ));
                  lickRatePostRewardTemp = noOfLicksPostReward/2; % Find LickRate for 2s of reward period
                  tempLickRateArray(tempPositionIndex) = lickRatePostRewardTemp;
              

                  if contains(data.(3)(currentRowNo), 'No Reward' )  % TrialType =1 for 'No Reward' trials, TrialType =2 for 'Reward' trials
                    trialType =1;
                  else 
                    trialType =2;
                  end     
                  
              end
           end
            
           tempPositionIndex =tempPositionIndex+1; % Next position
                   
        end
        
        currentRowNo =currentRowNo+1; % Next row
        
    end
    
    nextRowNo = currentRowNo-1; % As outer loop increases loop index, 1 is deducted here
     
 end 
 

 function [averageDataArray,data_SEM] = findMeanSEM(dataArray)
     % This functions returns Mean array and Standard Error of Mean    
 
     [noOfTrials,noOfPositions]=size(dataArray); 
   
     data_SEM = std(dataArray) /sqrt(noOfTrials); % Find Standard Error of Mean (SEM) 
     
     if noOfTrials <=1
       data_SEM= zeros(1, noOfPositions); % SEM is zero if only one trial is found
     end 
     
     if noOfTrials  > 1
        averageDataArray =  mean(dataArray);
     else 
        averageDataArray = dataArray;
     end 
 
 end 
