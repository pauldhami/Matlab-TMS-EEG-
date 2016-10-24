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

Subtract = true; %true, false, 'both'
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
% electrodes (columns)

SP_file_chan_count = pop_loadset('filename',...
        SP_list(1).name,...
        'filepath', pathin);

LICI_each_subject = zeros(size(SP_list,1),size(SP_file_chan_count.chanlocs(1,1:end),2));

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
    
    % Add individual subject info to structure
    TMSEEG(i).subjectID = SP_file.filename; 
    
    % Print row and file names to command window
    fprintf('\n For i (row) value %d, Single Pulse file is:\n\n %s \n\n and Paired Pulse file is:\n\n %s\n\n', i, SP_file.filename, PP_file.filename) 
    
    % Find the corresponding element index of desired times, to index EEG.data
    Time1_index = find(SP_file.times == Time1,1);
    Time2_index = find(SP_file.times == Time2,1); 

    % Preallocating vector to contain LICI value for each electrode of ON
    % SUBJECT
    
    %LICI holds the one inhibition value at each electrode
    LICI = zeros(1, numel(SP_data(:,1,1)));
    
    %SP_TEP and PP_TEP will hold mean TEP at each electrode
    SP_TEP = zeros(numel(SP_data(1,:,1)), numel(SP_data(:,1,1)));
    PP_LICI_TEP = zeros(numel(SP_data(1,:,1)), numel(SP_data(:,1,1)));

    % Looping over each electrode (h represents electrodes 1:60) 
    for h = 1:numel(SP_data(:,1,1))  
    
        % Mean across all epochs at electrode i for SP and PP
        SP = mean(SP_data(h,:,:), 3);
        PP = mean(PP_data(h,:,:), 3);
        
        %Transfering mean SP_TEP of each electrode transposed
        SP_TEP(:,h) = SP'; 
        
        % Apply subtraction method or not
        if Subtract == false
            
            if h == 1
                disp('No PP-SP subtraction');
            end
            
            PP_LICI_TEP(:,h) = PP';
            
            SP_abs = abs(SP);
            PP_abs = abs(PP);
                      
            SP_area = trapz(SP_abs(Time1_index:Time2_index));
            PP_area = trapz(PP_abs(Time1_index:Time2_index));
            
            LICI_elec = ((1 - (PP_area/SP_area))*100);
            LICI(h) = LICI_elec;
            LICI_each_subject(i,h) = LICI_elec;
            
            %LICI(1,i) = SP.chanlocs(1,i).labels;
            
        elseif Subtract == true
            
            if h == 1
                disp('PP-SP subtraction applied');
            end
            
            % Adjusting PP by lining SP TS (0-999ms) with PP CS (-100-899ms)
            % and subtracting
            
            PP_adj = PP;
            PP_adj(901:1900) = PP(901:1900) - SP(1001:2000);
            
            PP_LICI_TEP(:,h) = PP_adj';
            
            SP_abs = abs(SP);
            PP_adj_abs = abs(PP_adj);
            
            SP_area = trapz(SP_abs(Time1_index:Time2_index));
            PP_area = trapz(PP_adj_abs(Time1_index:Time2_index));
            
            LICI_elec = ((1 - (PP_area/SP_area))*100);
            LICI(h) = LICI_elec;
            LICI_each_subject(i,h) = LICI_elec;
            %LICI(1,i) = SP.chanlocs(1,i).labels;
            
        % elseif all(Subtract == 'both')
            
        else
            error('Improper argument for substract (true) or not (false). LICI not calculated');
            
        end   
    end
    
    % Allocating TEP info into participant 'i's structure
    
    TMSEEG(i).SP_TEP = SP_TEP;
    
    if Subtract == true
        TMSEEG(i).PP_TEP_subtraction = PP_LICI_TEP;
        TMSEEG(i).LICI_subtraction = LICI;
    elseif Subtract == false
        TMSEEG(i).PP_TEP_no_subtraction = PP_LICI_TEP;
        TMSEEG(i).LICI_no_subtraction = LICI;
    end
    
end

%%

% Averaging across subjects (average of each column in LICI_each_subject) 

LICI_subject_average = mean(LICI_each_subject, 1);

%%
% Topo plotting resulting AVERAGED ACROSS SUBJECTS LICI results 

figure
subplot(2,1,1)
topoplot(LICI_subject_average, SP_file.chanlocs); colorbar; caxis([-50 50])
title('LICI Inhibition (%) Topography Map', 'fontsize',16)

subplot(2,1,2)
hist(LICI_subject_average)
title('Count of LICI Inhibition Values (%) Across Electrodes', 'fontsize', 16) 

% FZ for subtraction and no subtraction for subject by subject with PP and
% SP in each figure

% Amplitude of N100 and correlation with LICI values and check how it
% changes with subtraction vs no subtraction 

