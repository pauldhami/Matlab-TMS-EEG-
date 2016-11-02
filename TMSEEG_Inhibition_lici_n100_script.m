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

% Try to commit


% October 26
% Source localization + some introduction of behavior (not even EEG, but 
% for example, NIH toolbox, GnG behaviour, WM etc.) 
% Maybe take average of sp and then subtract from pp trial by trial 
% MEAN TAKEN FIRST, THEN SUBTRACTION IS DONE?


% OCTOBER 24 :
% FZ for subtraction and no subtraction for subject by subject with PP and
% SP in each figure

% Amplitude of N100 and correlation with LICI values and check how it
% changes with subtraction vs no subtraction 



%%
Time1 = 50;
Time2 = 275;
pathin = '/Users/Prabh/Desktop/TMS_function_draft/files/RPFC_LICI/';

%%

% Creating structure array to hold individual subject info 

TMSEEG = struct('subjectID', {},...
    'SP_TEP', {},...
    'PP_TEP_subtraction', {},...
    'LICI_subtraction', {},...
    'PP_TEP_no_subtraction', {},...
    'LICI_no_subtraction', {},...
    'N100',{}); 


%%

% Pathin is location of files, SP and PP list contain file names of 
% SP and PP files found in pathin folder

SP_list = dir([pathin '*sp*']); 
 
PP_list = dir([pathin '*pp*']); 

if size(SP_list,1) ~= size(PP_list,1)
    error('Number of SP and PP files in folder not equal')
end

%%

% Getting channel info from single subject's dataset

SP_file_chan_count = pop_loadset('filename',...
        SP_list(1).name,...
        'filepath', pathin);
    
chaninfo = {SP_file_chan_count.chanlocs.labels};    
    
clear SP_file_chan_count

%%

% Preallocating vector to hold LICI data for each subject(rows) across
% each electrodes (columns)

LICI_each_subject_subtraction = zeros(size(SP_list,1),size(chaninfo,2)); 
LICI_each_subject_no_subtraction = zeros(size(SP_list,1),size(chaninfo,2));


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

% Averaging across subjects (average of each column in LICI_each_subject) 

LICI_subtraction_average = mean(LICI_each_subject_subtraction, 1);
LICI_no_subtraction_average = mean(LICI_each_subject_no_subtraction, 1); 

%% 
% N100 calculation for each subject, at each electrode 

N100_min_time = 90;
N100_max_time = 130; 

N100_each_subject_min = zeros(size(SP_list,1),size(chaninfo,2));

for i = 1:size(SP_list,1)
    
    SP_file = pop_loadset('filename',...
        SP_list(i).name,...
        'filepath', pathin);
    
    SP_data = SP_file.data;
    
    N100_subject_min = zeros(1, numel(SP_data(:,1,1)));
    
    for h = 1:numel(SP_data(:,1,1))
        
        SP = mean(SP_data(h,:,:), 3);
        N100_min_time_index = find(SP_file.times == N100_min_time);
        N100_max_time_index = find(SP_file.times == N100_max_time);
        N100_min_value = min(SP(N100_min_time_index:N100_max_time_index));
        N100_subject_min(h) = N100_min_value;
        N100_each_subject_min(i,h) = N100_min_value;
    end
    
    TMSEEG(i).N100 = N100_subject_min;
    
end

N100_min_average = mean(N100_each_subject_min, 1);

%%
%Clearing up workspace

% Keep SP_file due to it being referenced for channel info

clear LICI_elec ...
    PP_abs PP_adj PP_adj_abs PP_area PP_data PP_file ...
    PP_file_name PP_list SP SP_TEP TP_abs...
    SP_area SP_data SP_file_chan_count SP_list ...
    Time1_index Time2_index ans h i k ...
    LICI_no_subtraction LICI_subtraction ...
    PP_TEP_no_subtraction PP_TEP_subtraction SP_abs ...
    N100_min_time_index N100_max_time_index N100_min_value


%% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ~~~~~~~~ Graphs and Figures ~~~~~~~~~ %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Topo plotting resulting AVERAGED ACROSS SUBJECTS LICI results 

figure

%LICI No Subtraction
subplot(2,2,1)
topoplot(LICI_no_subtraction_average, SP_file.chanlocs); colorbar; caxis([-50 50])
title('LICI No Subtraction Inhibition (%)', 'fontsize',14)

subplot(2,2,3)
hist(LICI_no_subtraction_average)
title('Count of LICI Inhibition Values (%) Across Electrodes', 'fontsize', 14) 

%LICI Subtraction
subplot(2,2,2)
topoplot(LICI_subtraction_average, SP_file.chanlocs); colorbar; caxis([-50 50])
title('LICI Subtraction Applied Inhibition (%)', 'fontsize',14)

subplot(2,2,4)
hist(LICI_subtraction_average)
title('Count of LICI Inhibition Values (%) Across Electrodes', 'fontsize', 14)

%%

% Topo plotting resulting AVERAGED ACROSS SUBJECTS N100 results 

figure


subplot(2,1,1)
topoplot(N100_min_average, SP_file.chanlocs); colorbar; caxis([-6 6]) 
title('N100 Results Averaged Across Subjects', 'fontsize', 18) 

subplot(2,1,2)
hist(N100_min_average)



%%

% RECTIFIED (abs) waveform for LICI with/without subtraction, and SP

%Choose electrode

figure

electrode_to_plot = 'FCZ';

electrode_to_plot = find(strcmp(chaninfo,electrode_to_plot));

for i = 1:size(TMSEEG,2) 
    subplot(3,2,i)
    
    PP_TEP_no_sub = (abs(TMSEEG(i).PP_TEP_no_subtraction(:,electrode_to_plot))); 
    plot(PP_TEP_no_sub(1010:1300)) 
    hold on
    
    PP_TEP_sub = (abs(TMSEEG(i).PP_TEP_subtraction(:,electrode_to_plot)));
    plot(PP_TEP_sub(1010:1300), 'r') 
    
    hold on
    SP_TEP = (abs(TMSEEG(i).SP_TEP(:,electrode_to_plot))); 
    plot(SP_TEP(1010:1300), '--k')
    
    title(TMSEEG(i).subjectID(1:7), 'fontsize', 14)
end
legend('PP No Subtraction', 'PP Subtraction', 'SP'); 

%% 

%Each subjects paired topography maps for LICI with and without subtraction

for i = 1:6
    figure
    suptitle(TMSEEG(i).subjectID(1:7)) 
    
    subplot(2,2,1)
    topoplot(TMSEEG(i).LICI_subtraction, SP_file.chanlocs); colorbar
    title('LICI subtraction')
    caxis([-50 50])
    subplot(2,2,3)
    hist(TMSEEG(i).LICI_subtraction) 
    
    subplot(2,2,2)
    topoplot(TMSEEG(i).LICI_no_subtraction, SP_file.chanlocs); colorbar
    title('LICI no subtraction')
    caxis([-50 50])
    subplot(2,2,4)
    hist(TMSEEG(i).LICI_no_subtraction) 
end

%%

% TOPOPLOT OF ONLY SIGNIFICANT CORRELATIONS
elec_corr = zeros(1,60);

for elec = 1:60
    info_cor = zeros(2,size(TMSEEG,2));
   
    for i = 1:5
        info_cor(i,1) = TMSEEG(i).N100(elec); %N100 value
        info_cor(i,2) = TMSEEG(i).LICI_no_subtraction(elec);
    end
    
    [R, P] = corrcoef(info_cor(:,1), info_cor(:,2)); 
    if P(1,2) <= 0.05
        elec_corr(elec) = R(1,2);
    elseif P(1,2) > 0.05
        elec_corr(elec) = 0;
    end
end

figure
topoplot(elec_corr, SP_file.chanlocs)


%%

% CORRELATION TOPOPLOT

elec_corr = zeros(1,60);

% Correlation between N100 and LICI No Subtraction
for elec = 1:60
    info_cor = zeros(size(TMSEEG,2),2);
   
    for i = 1:size(TMSEEG,2)
        info_cor(i,1) = TMSEEG(i).N100(elec); %N100 value
        info_cor(i,2) = TMSEEG(i).LICI_no_subtraction(elec);
    end
    
    [R, P] = corrcoef(info_cor(:,1), info_cor(:,2)); 
   
    elec_corr(elec) = R(1,2);
end

figure
title('Correlation between N100 and LICI with no subtraction', 'fontsize',16)
colorbar; caxis([-1 1]) 
topoplot(elec_corr, SP_file.chanlocs)

% Correlation between N100 and LICI with Subtraction
for elec = 1:60
    info_cor = zeros(size(TMSEEG,2),2);
   
    for i = 1:size(TMSEEG,2)
        info_cor(i,1) = TMSEEG(i).N100(elec); %N100 value
        info_cor(i,2) = TMSEEG(i).LICI_subtraction(elec);
    end
    
    [R, P] = corrcoef(info_cor(:,1), info_cor(:,2)); 
   
    elec_corr(elec) = R(1,2);
end

figure
title('Correlation between N100 and LICI with subtraction', 'fontsize',16)
colorbar; caxis([-1 1]) 
topoplot(elec_corr, SP_file.chanlocs)

% Correlation between LICI with and without subtraction
for elec = 1:60
    info_cor = zeros(size(TMSEEG,2),2);
   
    for i = 1:size(TMSEEG,2)
        info_cor(i,1) = TMSEEG(i).LICI_no_subtraction(elec); %N100 value
        info_cor(i,2) = TMSEEG(i).LICI_subtraction(elec);
    end
    
    [R, P] = corrcoef(info_cor(:,1), info_cor(:,2)); 
   
    elec_corr(elec) = R(1,2);
end

figure
title('Correlation between LICI with and without subtraction', 'fontsize',16)
colorbar; caxis([-1 1]) 
topoplot(elec_corr, SP_file.chanlocs)

%%

% Correlation scatter plot at each electrode

electrode_to_plot = 'CZ';

electrode_to_plot = find(strcmp(chaninfo,electrode_to_plot));

my_corr = zeros(size(TMSEEG,2),2);

for i = 1:size(TMSEEG,2)
my_corr(i,1) = TMSEEG(i).LICI_subtraction(electrode_to_plot);
my_corr(i,2) = TMSEEG(i).N100(electrode_to_plot);
end

figure
scatter(my_corr(:,1), my_corr(:,2), 250)

% hold on 
% 
% for i = 1:size(TMSEEG,2)
% my_corr(i,1) = TMSEEG(i).LICI_no_subtraction(electrode_to_plot);
% my_corr(i,2) = TMSEEG(i).N100(electrode_to_plot);
% end
% 
% scatter(my_corr(:,1), my_corr(:,2))

%%


figure; 
pop_plottopo(TMSEEG(1).PP_TEP_no_subtraction, SP_file.chanlocs[1:60] , 0, 'ydir',1);


figure; pop_plottopo(TMSEEG, [1:60] , 'nothing', 0, 'ydir',1);

%%

average of sp and take from pp trial by trial 
