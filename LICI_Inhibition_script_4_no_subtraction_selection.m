%%%%%%%%%%%%%%%%%%%%%%%%%%
% LICI Inhibition Script %
%%%%%%%%%%%%%%%%%%%%%%%%%%

% Script that reads in SP and corresponding PP file and performs LICI 
% calculation at each electrode for each subject, then averages across
% subjects for a group average matrix/topo

% Indicate whether raw subtraction (true/false) and the 
% time range (Time1, Time2) that integration of rectified curve should 
% take place

% For now, make sure all types of files (location of stimulation and visit
% e.g. LPFC V14) are in a SEPERATE folder each, then change pathin to each
% folder and run script to create summary/figure of each factor level
% (location/visit)

% Still need to manually edit CS and TS timing for subtraction (line 120)

Time1 = 50;
Time2 = 275;
pathin = '/Users/Prabh/Desktop/TMS_function_draft/files/LPFC_LICI/';

%%

% Pathin is location of files, SP and PP list contain file names of 
% SP and PP files found in pathin folder

SP_list = dir([pathin '*sp*']); 
 
PP_list = dir([pathin '*pp*']); 

if size(SP_list,1) ~= size(PP_list,1)
    error('Number of SP and PP files in folder not equal')
end

%%

% Preallocating vector to hold LICI data for each subject(rows) across
% each electrodes (columns)

SP_file_chan_count = pop_loadset('filename',...
        SP_list(1).name,...
        'filepath', pathin);

LICI_each_subject_subtraction = zeros(size(SP_list,1),size(SP_file_chan_count.chanlocs(1,1:end),2));
LICI_each_subject_no_subtraction = zeros(size(SP_list,1),size(SP_file_chan_count.chanlocs(1,1:end),2));

%%

% Creating structure array to hold individual subject info 

TMSEEG = struct('subjectID', {},...
    'SP_TEP', {},...
    'PP_TEP_subtraction', {},...
    'LICI_subtraction', {},...
    'PP_TEP_no_subtraction', {},...
    'LICI_no_subtraction', {}); 


%%

% Looping over each pair of SP and PP files (looping through each subject)

for i = 1:size(SP_list,1)
    
    % Selecting SP file
    SP_file = pop_loadset('filename',...
        SP_list(i).name,...
        'filepath', pathin);
    
    SP_data = SP_file.data;
    
    % Selecting PP file that is matched to SP file chosen above
    PP_file_name = [];
    for k = 1:size(PP_list)
        if (all(PP_list(k).name(1:7) == SP_file.filename(1:7)) && ~(strcmp(PP_list(k).name, SP_file.filename))) 
            PP_file_name = PP_list(k).name;
        end  
    end
    
    PP_file = pop_loadset('filename',...
        PP_file_name,...
        'filepath', pathin);
    
    PP_data = PP_file.data;
    
    % Errors thrown if length of SP and PP data are different; different
    % electrode count
    if SP_file.pnts ~= PP_file.pnts 
    error('SP and PP are not same length')
    end

    if all(SP_file.times ~= PP_file.times)
    error('SP and PP have different time points')
    end
    
    if any([PP_file.chanlocs(1,:).labels] ~= [SP_file.chanlocs(1,:).labels])
    error('PP and SP not same elec info')
    end
    
    % Add individual subject ID to structure
    TMSEEG(i).subjectID = SP_file.filename; 
    
    % Print row and file names to command window
    fprintf('\n For i (row) value %d, Single Pulse file is:\n\n %s \n\n and Paired Pulse file is:\n\n %s\n\n', i, SP_file.filename, PP_file.filename) 
    
    % Find the corresponding element index of desired times, to index EEG.data
    Time1_index = find(SP_file.times == Time1,1);
    Time2_index = find(SP_file.times == Time2,1); 

    % Preallocating vector to contain LICI value for each electrode of ONE
    % SUBJECT
    
    %LICI holds the one inhibition value at each electrode
    LICI_subtraction = zeros(1, numel(SP_data(:,1,1)));
    LICI_no_subtraction = zeros(1, numel(SP_data(:,1,1)));
    
    %SP_TEP and PP_TEP will hold mean TEP at each electrode
    SP_TEP = zeros(numel(SP_data(1,:,1)), numel(SP_data(:,1,1)));
    PP_TEP_subtraction = zeros(numel(SP_data(1,:,1)), numel(SP_data(:,1,1)));
    PP_TEP_no_subtraction = zeros(numel(SP_data(1,:,1)), numel(SP_data(:,1,1)));

    % Looping over each electrode (h represents electrodes 1:60) 
    for h = 1:numel(SP_data(:,1,1))  
    
        % Mean across all epochs at electrode h for SP and PP
        SP = mean(SP_data(h,:,:), 3);
        PP = mean(PP_data(h,:,:), 3);
        
        %Transfering mean SP_TEP of each electrode transposed
        SP_TEP(:,h) = SP'; 
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % No subtraction applied for LICI %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
        PP_TEP_no_subtraction(:,h) = PP';
            
        SP_abs = abs(SP);
        PP_abs = abs(PP);
                      
        SP_area = trapz(SP_abs(Time1_index:Time2_index));
        PP_area = trapz(PP_abs(Time1_index:Time2_index));
            
        LICI_elec = ((1 - (PP_area/SP_area))*100);
        LICI_no_subtraction(h) = LICI_elec;
        LICI_each_subject_no_subtraction(i,h) = LICI_elec;
            
            %LICI(1,i) = SP.chanlocs(1,i).labels;
            
        %%%%%%%%%%%%%%%%%%%%%%%    
        % Subtraction applied %
        %%%%%%%%%%%%%%%%%%%%%%%
        
        % Adjusting PP by lining SP TS (0-999ms) with PP CS (-100-899ms)
        % and subtracting
            
        PP_adj = PP;
        PP_adj(901:1900) = PP(901:1900) - SP(1001:2000);
            
        PP_TEP_subtraction(:,h) = PP_adj';
            
        SP_abs = abs(SP);
        PP_adj_abs = abs(PP_adj);
            
        SP_area = trapz(SP_abs(Time1_index:Time2_index));
        PP_area = trapz(PP_adj_abs(Time1_index:Time2_index));
            
        LICI_elec = ((1 - (PP_area/SP_area))*100);
        LICI_subtraction(h) = LICI_elec;
        LICI_each_subject_subtraction(i,h) = LICI_elec;
        %LICI(1,i) = SP.chanlocs(1,i).labels;
            
    end
    
    % Allocating TEP info into participant 'i's structure
    
    TMSEEG(i).SP_TEP = SP_TEP;
    
    
    TMSEEG(i).PP_TEP_subtraction = PP_TEP_subtraction;
    TMSEEG(i).LICI_subtraction = LICI_subtraction;
  
    TMSEEG(i).PP_TEP_no_subtraction = PP_TEP_no_subtraction;
    TMSEEG(i).LICI_no_subtraction = LICI_no_subtraction;
  
    
end

%%

%Clearing up workspace

% Keep SP_file due to it being referenced for channel info

clear LICI_elec PP_abs PP_adj PP_adj_abs PP_area PP_data PP_file ...
    PP_file_name PP_list SP SP_TEP TP_abs SP_area SP_data ...
    SP_file_chan_count SP_list Time1_index Time2_index ans h i k ...
    LICI_no_subtraction LICI_subtraction ...
    PP_TEP_no_subtraction PP_TEP_subtraction SP_abs
    

%%

% Averaging across subjects (average of each column in LICI_each_subject) 

LICI_subtraction_average = mean(LICI_each_subject_subtraction, 1);
LICI_no_subtraction_average = mean(LICI_each_subject_no_subtraction, 1); 

%%
% Topo plotting resulting AVERAGED ACROSS SUBJECTS LICI results 


avg = LICI_subtraction_average; 

figure
subplot(2,1,1)
topoplot(avg, SP_file.chanlocs); colorbar; caxis([-50 50])
title('LICI Inhibition (%) Topography Map', 'fontsize',16)

subplot(2,1,2)
hist(avg)
title('Count of LICI Inhibition Values (%) Across Electrodes', 'fontsize', 16) 

% FZ for subtraction and no subtraction for subject by subject with PP and
% SP in each figure

% Amplitude of N100 and correlation with LICI values and check how it
% changes with subtraction vs no subtraction 




%%

figure

for i = 1:6
subplot(3,2,i)
PP_TEP_no_sub = (abs(TMSEEG(i).PP_TEP_no_subtraction(:,10))); 
plot(PP_TEP_no_sub([1010:1275])) 
hold on
PP_TEP_sub = (abs(TMSEEG(i).PP_TEP_subtraction(:,10)));
plot(PP_TEP_sub([1010:1275]), 'r')  
title(TMSEEG(i).subjectID(1:7))
end

%%



for i = 1:6
    figure 
    subplot(1,2,1)
    topoplot(TMSEEG(i).LICI_subtraction, SP_file.chanlocs)
    title('LICI subtraction')
    caxis([-50 50])
    subplot(1,2,2)
    topoplot(TMSEEG(i).LICI_no_subtraction, SP_file.chanlocs)
    title('LICI no subtraction')
    caxis([-50 50])
end

