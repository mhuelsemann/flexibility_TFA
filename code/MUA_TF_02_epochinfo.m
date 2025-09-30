%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Title:         mass-univariate analysis with cluster-based permutation  %
%                tests of time-frequency data                             %
% Subtitle:      extracting epoch infos for task and                      %
%                find participants with unsufficient task performance and %
%                summarise preprocessing stats                            %
% Author:        Mareike J. Hülsemann, University of Mainz                %
% Date:          11/2023                                                  %
% Last Modified: 09/2025                                                  %
% Note:          Exp23, Switching Tasks                                   %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% BSD 3-Clause License
% 
% Copyright (c) 2025, Mareike Huelsemann
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
% 
% 1. Redistributions of source code must retain the above copyright notice, this
%    list of conditions and the following disclaimer.
% 
% 2. Redistributions in binary form must reproduce the above copyright notice,
%    this list of conditions and the following disclaimer in the documentation
%    and/or other materials provided with the distribution.
% 
% 3. Neither the name of the copyright holder nor the names of its
%    contributors may be used to endorse or promote products derived from
%    this software without specific prior written permission.
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
% FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
% DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
% SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
% CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
% OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
% OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.



%% ---------- Set Up ------------------------------------------------------

% define experiment and task specific settings
task_name     = 'oddeven'; % 'oddeven', 'globallocal', 'numberletter'
min_correct   = .7;

% define PATHs
[filepath,~,~] = fileparts(matlab.desktop.editor.getActiveFilename);
parts          = regexp(filepath,filesep,'split');
PATH_MAIN      = fullfile(parts{1:end-1});

PATH_EEGLAB    = fullfile(PATH_MAIN, 'eeglab2023.0');

PATH_MAIN   = fullfile(PATH_MAIN, 'derivatives', task_name);
PATH_IN     = fullfile(PATH_MAIN, 'preprocessed');
PATH_STATS  = fullfile(PATH_MAIN, 'stats');
PATH_OUT    = fullfile(PATH_MAIN, 'stats', 'summary');

% load list of continuous preprocessed data sets
allSets = dir(fullfile(PATH_IN, '*.set'));


%% ---------- Preprocessing Stats ------------------------------------------

% initialize preprocessing stats
artefact_varnames = {...
    'id', ...
    'task', ...
    'nb_chan_raw', ...
    'pnts_raw', ...
    'nb_events_raw', ...
    'min_num_pnts_ICA_20', ...
    'min_num_pnts_ICA_30', ...
    'pnts_task', ...
    'nb_events_task', ...
    'nb_trials_task', ...
    'srate_raw', ...
    'task_length_min', ...
    'srate_for_ICA', ...
    'nb_chan_cleaned', ...
    'nb_trials', ...
    'pnts_trial', ...
    'nb_trials_clean', ...
    'nb_trials_rejected', ...
    'prct_trials_kept', ...
    'prct_trials_rejected', ...
    'nb_trials_rejected_total', ...
    'nb_trials_rejected_thre', ...
    'nb_trials_rejected_prob', ...
    'nb_trials_rejected_kurt', ...
    'nb_trials_rejected_unique_thre', ...
    'nb_trials_rejected_unique_prob', ...
    'nb_trials_rejected_unique_kurt', ...
    'nb_trials_rejected_inTP_notK', ...
    'nb_trials_rejected_inTK_notP', ...
    'nb_trials_rejected_inPK_notT', ...
    'nb_trials_rejected_inTPK', ...
    'seed_ica', ...
    'time_min_runica', ...
    'nb_IC_brain', ...
    'nb_IC_muscle', ...
    'nb_IC_eye', ...
    'nb_IC_heart', ...
    'nb_IC_linenoise', ...
    'nb_IC_channoise', ...
    'nb_IC_other', ...
    'median_prct_IC_brain', ...
    'median_prct_IC_muscle', ...
    'median_prct_IC_eye', ...
    'median_prct_IC_heart', ...
    'median_prct_IC_linenoise', ...
    'median_prct_IC_channoise', ...
    'median_prct_IC_other', ...
    'nb_rejIC_eye', ...
    'nb_rejIC_muscle', ...
    'nb_rejIC_heart', ...
    'nb_rejIC_total'};
artefact_overview  = nan(numel(allSets), numel(artefact_varnames));
channel_overview   = nan(numel(allSets), 32);

% summarise individual stats
switch task_name
    case 'oddeven',      letters2remove = 35;
    case 'globallocal',  letters2remove = 39;
    case 'numberletter', letters2remove = 40;
    otherwise
        warning('Typo in variable [task_name].')
end
for n = 1:numel(allSets)
   % artefacts
   tmp = load(fullfile(PATH_STATS, [allSets(n).name(1:end-letters2remove) 'preprocessing_artefact_overview.csv']));
   artefact_overview(n,:) = tmp;    
   
   % channels
   tmp = load(fullfile(PATH_STATS, [allSets(n).name(1:end-letters2remove) 'preprocessing_artefact_chans.csv']));
   % tmp(isnan(tmp)) = [];
   chans = zeros(32,1);
   chans(setdiff(1:32, tmp)) = 1;   
   channel_overview(n,:) = chans;
end

% save summarised stats
artefact_overview = array2table(artefact_overview);
artefact_overview.Properties.VariableNames = artefact_varnames;
artefact_overview = addvars(artefact_overview, char({allSets.name}'), 'Before', 'id');
writetable(artefact_overview, fullfile(PATH_OUT, [allSets(1).name(1:6) task_name '_preprocessing_artefact_overview_all.csv']))

channel_overview = array2table(channel_overview);
channel_overview.Properties.VariableNames = string(num2cell(1:32));
channel_overview = addvars(channel_overview, char({allSets.name}'), 'Before', '1');
writetable(channel_overview, fullfile(PATH_OUT, [allSets(1).name(1:6) task_name '_preprocessing_artefact_chans_all.csv']))


%% ---------- start EEGLAB ------------------------------------------------
addpath(PATH_EEGLAB);
eeglab; close(gcf);
pop_editoptions('option_single', 0, ...
    'option_savetwofiles', 0, ... % new default: FALSE
    'option_computeica',   1, ... % one only has EEG.icaact if TRUE
    'option_checkversion', 0);


%% ---------- Add epoch infos ---------------------------------------------

% initialise
epoch_info_varnames = {'Trial', 'Task', 'Switch', 'RT', 'Accuracy', ...
    'miss', 'z_log_RT', 'outlier', 'rel_odd_even', 'EEGreject'};
trialscounts = artefact_overview.nb_trials_task;
task_performance = table(char({allSets.name}'), nan(numel(allSets),1));
task_performance.Properties.VariableNames = {'Set', 'Acc'};
exclude = [];
    
% open file to save exlusion due to insufficient task performance
fileID = fopen(fullfile(PATH_OUT, ...
            [allSets(1).name(1:6) 'exclude_task_performance.txt']),'w');


for n=1:length(allSets)
    disp(['**************** Started Set ' num2str(n) ' of ' num2str(length(allSets)) ' ****************'])
    
    % initialise
    trialcount = trialscounts(n);
    switch task_name
        case 'oddeven'
            if n==2, trialcount = trialcount - 1; end % remove last incomplete trial
        case 'numberletter'
            if n==110, trialcount = trialcount - 1; end
    end
    epoch_info = nan(trialcount, numel(epoch_info_varnames));
    
    % load data set
    EEG = pop_loadset('filename',allSets(n).name,'filepath',PATH_IN);
    
    % prepare event data
    stim    = {EEG.event.type};
    latency = [EEG.event.latency];
    
    if strcmp(task_name,'numberletter') && n==110
        event1 = 16;
    else
        event1 = find(strcmp(stim, 'S101')) + 1;
    end

    if (strcmp(task_name,'oddeven') && n==2) || ...
            (strcmp(task_name,'numberletter') && (n==56 || n==66 || n==72))
        event2 = event1 + (trialcount)*5 - 1;
    else
        event2 = find(strcmp(stim, 'S201')) - 1;
    end
    
    stim    = stim(event1:event2);
    latency = latency(event1:event2);
    
    switch task_name
        case 'oddeven'
            if n == 9
                list1 = strcmp(stim, 'S255');
                list1 = [list1(end) list1(1:end-1)];
                list2 = strcmp(stim, 'S250');
                marker2del = find(list1 & list2);
                stim(marker2del) = [];
                latency(marker2del) = [];
            end
    end

    stim    = reshape(stim,5,trialcount);
    latency = reshape(latency,5,trialcount);    
    
    % sanity check
    if sum(matches(stim(1,:), {'S 31', 'S 32', 'S 33', 'S 34'})) ~= trialcount, warning(['Set ' num2str(n) ': ' allSets(n).name ' Error for fix cross marker (S 3*)']), end
    if sum(matches(stim(2,:), {'S 41', 'S 42', 'S 43', 'S 44'})) ~= trialcount, warning(['Set ' num2str(n) ': ' allSets(n).name ' Error for ISI marker (S 4*)']), end
    if sum(matches(stim(3,:), {'S 51', 'S 52', 'S 53', 'S 54'})) ~= trialcount, warning(['Set ' num2str(n) ': ' allSets(n).name ' Error for target marker (S 5*)']), end
    if sum(matches(stim(4,:), {'S150', 'S160', 'S250', 'S251', 'S255', 'S222'})) ~= trialcount, warning(['Set ' num2str(n) ': ' allSets(n).name ' Error for response marker (S1**/S2**)']), end
    if sum(matches(stim(5,:), {'S 90'})) ~= trialcount, warning(['Set ' num2str(n) ': ' allSets(n).name ' Error for ITI marker (S 90)']), end
    
    % calculate RT
    rt = latency(4,:) - latency(3,:);
    
    % define conditions
    col=    1; epoch_info(:,col) = 1:size(epoch_info,1); % #rowtrial
    col=col+1; % task {1=LessMore; 2=OddEven}
    epoch_info(ismember(stim(1,:),{'S 31' 'S 32'}), col) = 1;
    epoch_info(ismember(stim(1,:),{'S 33' 'S 34'}), col) = 2;
    col=col+1; % switch {1=Switch; 0=Repeat}
    epoch_info(ismember(stim(1,:),{'S 31' 'S 33'}),col) = 0;
    epoch_info(ismember(stim(1,:),{'S 32' 'S 34'}),col) = 1;
    col=col+1; % RT
    epoch_info(:,col) = rt; 
    col=col+1; % accuracy {1=Correct; 0=Wrong}
    epoch_info(ismember(stim(4,:),{'S150' 'S160'}),col) = 1;
    epoch_info(ismember(stim(4,:),{'S250' 'S255' 'S251' 'S222'}),col) = 0;
    col=col+1;
    epoch_info(:,col) = 0; % response given
    epoch_info(ismember(stim(4,:),{'S222'}),col) = 1; % miss
    
    % outlier detection (+/- 3 SD) for z-scaled log(RT), condition specific
    col=col+1;
    rows = epoch_info(:,2) == 1 & epoch_info(:,3) == 0 & epoch_info(:,5) == 1; % LessMore-Repeat (only correct trials)
    epoch_info(rows,col) = zscore(log(rt(rows))); 
    rows = epoch_info(:,2) == 1 & epoch_info(:,3) == 1 & epoch_info(:,5) == 1; % LessMore-Switch (only correct trials)
    epoch_info(rows,col) = zscore(log(rt(rows))); 
    rows = epoch_info(:,2) == 2 & epoch_info(:,3) == 0 & epoch_info(:,5) == 1; % OddEven-Repeat (only correct trials)
    epoch_info(rows,col) = zscore(log(rt(rows))); 
    rows = epoch_info(:,2) == 2 & epoch_info(:,3) == 1 & epoch_info(:,5) == 1; % OddEven-Switch (only correct trials)
    epoch_info(rows,col) = zscore(log(rt(rows))); 
    col=col+1; % outlier
    epoch_info(:,col) = 0; % initial non-outlier
    epoch_info(:,col) = abs(epoch_info(:,col-1)) > 3; % outlier
    
    % odd-even split for reliability
    vec = repmat([0 1],1,48);
    col=col+1;
    rows = epoch_info(:,2) == 1 & epoch_info(:,3) == 0 & epoch_info(:,5) == 1 & epoch_info(:,8) == 0; % LessMore-Repeat (correct, without outliers)
    epoch_info(rows,col) = vec(1:sum(rows));
    rows = epoch_info(:,2) == 1 & epoch_info(:,3) == 1 & epoch_info(:,5) == 1 & epoch_info(:,8) == 0; % LessMore-Switch (correct, without outliers)
    epoch_info(rows,col) = vec(1:sum(rows));
    rows = epoch_info(:,2) == 2 & epoch_info(:,3) == 0 & epoch_info(:,5) == 1 & epoch_info(:,8) == 0; % OddEven-Repeat (correct, without outliers)
    epoch_info(rows,col) = vec(1:sum(rows));
    rows = epoch_info(:,2) == 2 & epoch_info(:,3) == 1 & epoch_info(:,5) == 1 & epoch_info(:,8) == 0; % OddEven-Switch (correct, without outliers)
    epoch_info(rows,col) = vec(1:sum(rows));
    
    % sanity check
    if sum(isnan(epoch_info(:,2))), warning(['Set ' num2str(n) ': ' allSets(n).name ' Error for task']), end
    if sum(isnan(epoch_info(:,3))), warning(['Set ' num2str(n) ': ' allSets(n).name ' Error for switch']), end
    if sum(isnan(epoch_info(:,4))), warning(['Set ' num2str(n) ': ' allSets(n).name ' Error for rt']), end
    if sum(isnan(epoch_info(:,5))), warning(['Set ' num2str(n) ': ' allSets(n).name ' Error for acc']), end
    if sum(isnan(epoch_info(:,6))), warning(['Set ' num2str(n) ': ' allSets(n).name ' Error for miss']), end
    if sum(isnan(epoch_info(:,8))), warning(['Set ' num2str(n) ': ' allSets(n).name ' Error for outlier']), end
    
    % add epoch info to EEG data set
    EEG.etc.epoch_info = epoch_info;
    EEG.etc.epoch_info_varnames = epoch_info_varnames;
    
    % save EEG dataset
    EEG = pop_saveset(EEG, 'filename', allSets(n).name, 'filepath', PATH_IN);
    
    % calculate task performance
    accuracy = mean(EEG.etc.epoch_info(:,5), 'omitna');
    task_performance{n,2} = accuracy;
    
    if accuracy < min_correct
        fprintf(fileID, ...
            'Set %d (%s) -> Exclude as accuracy is %1.2f.\n', ...
            n, allSets(n).name, accuracy);
        exclude = [exclude n];
    end    
end
fclose(fileID);

% delete participants with insufficients task performance
artefact_overview(exclude,:) = [];
channel_overview(exclude,:)  = [];
task_performance(exclude,:)  = [];

% save individual statistics
writetable(artefact_overview,fullfile(PATH_OUT, [allSets(1).name(1:6) task_name '_preprocessing_artefact_overview_final.csv']))
writetable(channel_overview, fullfile(PATH_OUT, [allSets(1).name(1:6) task_name '_preprocessing_artefact_chans_final.csv']))
writetable(task_performance, fullfile(PATH_OUT, [allSets(1).name(1:6) task_name '_task_performance_final.csv']))

% calculate summarised statistics
artefact_summary  = nan(5, numel(artefact_varnames));

artefact_summary(1,3:end) = min(artefact_overview{:,4:end}, [], 'omitnan');
artefact_summary(2,3:end) = max(artefact_overview{:,4:end}, [], 'omitnan');
artefact_summary(3,3:end) = median(artefact_overview{:,4:end}, 'omitnan');
artefact_summary(4,3:end) = mean(artefact_overview{:,4:end}, 'omitnan');
artefact_summary(5,3:end) = std(artefact_overview{:,4:end}, 'omitnan');

artefact_summary = array2table(artefact_summary);
artefact_summary.Properties.VariableNames = artefact_varnames;
artefact_summary.id = char({'min'; 'max'; 'median'; 'mean'; 'sd'});
artefact_summary.task = repmat(task_name, 5, 1);

task_performance_summary = table(char({'min'; 'max'; 'median'; 'mean'; 'sd'}), ...
    [min(task_performance{:,2}, [], 'omitnan'); ...
    max(task_performance{:,2}, [], 'omitnan'); ...
    median(task_performance{:,2}, 'omitnan'); ...
    mean(task_performance{:,2}, 'omitnan'); ...
    std(task_performance{:,2}, 'omitnan')]);
task_performance_summary.Properties.VariableNames = {'id', 'acc'};    

% save summarised statistics
writetable(artefact_summary,fullfile(PATH_OUT, [allSets(1).name(1:6) task_name '_preprocessing_artefact_overview_final_summary.csv']))
csvwrite(fullfile(PATH_OUT, ...
    [allSets(1).name(1:6) task_name '_preprocessing_artefact_chans_final_summary.csv']), ...
    sum(channel_overview{:,2:end}))
writetable(task_performance_summary,fullfile(PATH_OUT, [allSets(1).name(1:6) task_name '_task_performance_final_summary.csv']))

