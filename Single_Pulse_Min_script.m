%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TMS-EEG Single Pulse (N100 - MIN) Script %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Pathin of SP data
pathin = '/Users/Prabh/Desktop/TMS_function_draft/files/LPFC_LICI/';
SP_list = dir([pathin '*sp*']);

% Loading single SP to get channel information
SP_file_chan_count = pop_loadset('filename',...
        SP_list(1).name,...
        'filepath', pathin);
    
% Structure array to hold individual subject info
TMSEEG = struct('N100', {}); 

% Setting N100 epoch     
N100_min_time = 90;
N100_max_time = 130; 

% Preallocating vector to hold individual N100 data (rows = subject,
% columns = electrodes)
N100_each_subject_min = zeros(size(SP_list,1),size(SP_file_chan_count.chanlocs(1,1:end),2));

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

% Collapsing N100_each_subject_min to get average N100 
N100_min_average = mean(N100_each_subject_min, 1);