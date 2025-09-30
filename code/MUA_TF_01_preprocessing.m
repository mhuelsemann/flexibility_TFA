%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Title:         mass-univariate analysis with cluster-based permutation  %
%                tests of time-frequency data                             %
% Subtitle:      preprocessing                                            %
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
switch task_name
    case 'oddeven'
        marker_start  = 'S  6';
        marker_end    = 'S  7';
        marker_next   = 'S115';
    case 'globallocal'        
        marker_start  = 'S 14';
        marker_end    = 'S 15';
        marker_next   = 'S114';
    case 'numberletter'
        marker_start  = 'S 12';
        marker_end    = 'S 13';
        marker_before = 'S116';
        marker_next   = 'S119';      
    otherwise
        warning('Typo in variable [task_name].')
end
marker_start_num = str2double(marker_start(3:end));
event_indices = {'S 31' 'S 32' 'S 33' 'S 34'};
online_ref    = struct('labels',{'Cz'}, ...
    'sph_radius', {1}, 'sph_theta', {0}, 'sph_phi', {0}, ...
    'theta',{0}, 'radius',{0}, ...
    'X',{0},'Y',{0},'Z',{1}, ...
    'type',{''},'ref',{''},'urchan',{[]}); %,'datachan',{0});

% define PATHs
[filepath,~,~] = fileparts(matlab.desktop.editor.getActiveFilename);
parts          = regexp(filepath,filesep,'split');
PATH_MAIN      = fullfile(parts{1:end-1});

PATH_EEGLAB     = fullfile(PATH_MAIN, 'eeglab2023.0');

PATH_IN         = fullfile(PATH_MAIN, 'data');
PATH_OUT        = fullfile(PATH_MAIN, 'derivatives', task_name);
PATH_OUT_STATS  = fullfile(PATH_OUT, 'stats');
PATH_OUT_LN     = fullfile(PATH_OUT_STATS, 'LN');
PATH_OUT_ICA    = fullfile(PATH_OUT_STATS, 'ICA');
PATH_OUT_ICL    = fullfile(PATH_OUT_STATS, 'ICL');
PATH_OUT_PREPRO = fullfile(PATH_OUT, 'preprocessed');

% define FILE names
FILE_OUT            = ['_' task_name '_preprocessed_continuous'];
FILE_STATS_OVERVIEW = '_preprocessing_artefact_overview';
FILE_STATS_EPOCHS   = '_preprocessing_artefact_epochs';
FILE_STATS_CHANS    = '_preprocessing_artefact_chans';

% load list of raw EEG data
allVhdr = dir(fullfile(PATH_IN, '*.vhdr'));
% refine subject list
switch task_name
    case 'oddeven'
        files2exclude = {'Exp23_0043.vhdr', 'Exp23_0113.vhdr', ...
            'Exp23_0113_1.vhdr', 'Exp23_0119.vhdr'};
    case 'globallocal'
        files2exclude = {'Exp23_0043_2.vhdr', 'Exp23_0113_1.vhdr', ...
            'Exp23_01191.vhdr'};
    case 'numberletter'
        files2exclude = {'Exp23_0006_2_2.vhdr', 'Exp23_0006_2_3.vhdr', ...
            'Exp23_0050_2.vhdr', 'Exp23_0050_22.vhdr', ...
            'Exp23_0050_222.vhdr', 'Exp23_0052_2.vhdr', ...
            'Exp23_0063_22.vhdr', 'Exp23_0063_23.vhdr', ...
            'Exp23_0073_2_2.vhdr', 'Exp23_0080_2_2.vhdr', ...
            'Exp23_0119_2_2.vhdr', 'Exp23_0135_2_Posner.vhdr', ...
            'Exp23_0138_21.vhdr'};
    otherwise
        warning('Typo in variable [task_name].')
end
idx2exclude = find(ismember({allVhdr.name}, files2exclude));
allVhdr(idx2exclude) = [];

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
artefact_epochs_varnames = {'trials_rejected_unique_thre', ...
    'trials_rejected_unique_prob', ...
    'trials_rejected_unique_kurt', ...
    'trials_rejected_inTP_notK', ...
    'trials_rejected_inTK_notP', ...
    'trials_rejected_inPK_notT', ...
    'trials_rejected_inTPK', ...
    'trials_rejected_all'};


%% ---------- start EEGLAB ------------------------------------------------
addpath(PATH_EEGLAB);
eeglab; close(gcf);
pop_editoptions('option_single', 0, ...
    'option_savetwofiles', 0, ... % new default: FALSE
    'option_computeica',   1, ... % one only has EEG.icaact if TRUE
    'option_checkversion', 0);


%% ---------- Loop Over Subjects ------------------------------------------
for n=1:length(allVhdr)
    disp(['**************** Started Set ' num2str(n) ' of ' num2str(length(allVhdr)) '; ' allVhdr(n).name '; ' datestr(datetime()) ' ****************'])
    
    chans              = cell(32,1); % save channels after ASR procedure
    artefact_overview  = nan(1,  numel(artefact_varnames));
    artefact_epochs    = cell(1, numel(artefact_epochs_varnames));


    % ---------- Loading data set -----------------------------------------
    % load data set
    EEG = pop_loadbv(PATH_IN, allVhdr(n).name, [], []);
    
    % remove unused channels (e. g. EOG, ECG)
    EEG = pop_select( EEG, 'nochannel',{'vEOG_o' 'vEOG_u'});
    
    % save urchanlocs structure with online reference
    orig_chanlocs = EEG.chanlocs;
    orig_chanlocs(EEG.nbchan+1) = online_ref;
    EEG.etc.orig_chanlocs = orig_chanlocs;
    
    % document data set statistics
    col=    1; artefact_overview(1,col) = str2double(allVhdr(n).name(7:end-5));
    col=col+1; artefact_overview(1,col) = marker_start_num;
    col=col+1; artefact_overview(1,col) = EEG.nbchan;
    col=col+1; artefact_overview(1,col) = EEG.pnts;
    col=col+1; artefact_overview(1,col) = numel(EEG.event);
    
    % segment task with sufficient buffer time for filtering; e.g. +/- 1min
    buffer    = 30; % seconds    
    try    TaskStart = EEG.event(find(strcmp({EEG.event.type}, marker_start))).latency;
    catch
        try
            % use "previous task start"-marker instead of "current task start"-marker
            alternative_marker = find(strcmp({EEG.event.type}, marker_before));
            TaskStart   = EEG.event(alternative_marker(1)).latency;
        catch
            TaskStart = 1;
        end
    end  
    try    TaskEnd   = EEG.event(find(strcmp({EEG.event.type}, marker_end))).latency;
    catch
        try 
            % use "next task start"-marker instead of "current task end"-marker
            alternative_marker = find(strcmp({EEG.event.type}, marker_next));
            TaskEnd   = EEG.event(alternative_marker(1)).latency;
        catch
            TaskEnd   = EEG.pnts;
        end
    end
    TaskStart = max(TaskStart - buffer*EEG.srate, 1);
    TaskEnd   = min(TaskEnd   + buffer*EEG.srate, EEG.pnts);
    if TaskStart == 1 || TaskEnd == EEG.pnts
        fileID = fopen(fullfile(PATH_OUT_STATS, ...
            [allVhdr(n).name(1:end-5) '_warning.txt']),'w');
        fprintf(fileID,'Set %d (%s) -> There might not be enough buffer time for filtering or missing markers.',n, allVhdr(n).name(1:end-5));
        fprintf(fileID,'This warning was produced by >> MUA_TF_01_preprocessing.m');
        fclose(fileID);
    end
    EEG = pop_select(EEG,'point',[TaskStart TaskEnd]);

    % document data set statistics
    col=col+1; artefact_overview(1,col) = EEG.nbchan^2 * 20;
    col=col+1; artefact_overview(1,col) = EEG.nbchan^2 * 30;
    col=col+1; artefact_overview(1,col) = EEG.pnts;
    col=col+1; artefact_overview(1,col) = numel(EEG.event);
    col=col+1; artefact_overview(1,col) = sum(ismember({EEG.event.type },event_indices));    
    col=col+1; artefact_overview(1,col) = EEG.srate;   
    col=col+1; artefact_overview(1,col) = EEG.pnts/EEG.srate/60;   
    
    
    % ---------- ICA ------------------------------------------------------
    
    % create data set for ICA
    ICA = EEG;
    
    % document data set statistics
    col=col+1; artefact_overview(1,col) = ICA.srate;
        
    % filter data (1 Hz high pass for ICA)
    ICA = pop_eegfiltnew(ICA, 'locutoff', 1); % optional: 'hicutoff'
            
    % remove noisy channels
    ICA = clean_artifacts(ICA, ...
        'FlatlineCriterion','off', ...
        'ChannelCriterion',0.8, ...
        'LineNoiseCriterion','off', ...
        'Highpass','off', ...
        'BurstCriterion','off', ...
        'WindowCriterion','off', ...
        'BurstRejection','off', ...
        'Distance','Euclidian');
    chans(1:ICA.nbchan,1) = {ICA.chanlocs.labels};
    chans = cellfun(@str2double,chans);
    csvwrite(fullfile(PATH_OUT_STATS, ...
        [allVhdr(n).name(1:end-5), FILE_STATS_CHANS, '.csv']), chans)
    % document data set statistics
    col=col+1; artefact_overview(1,col) = ICA.nbchan;
    
    % remove line noise
    figure('Name', allVhdr(n).name(1:end-5),'NumberTitle','off', ...
        'Units', 'centimeters', 'Position', [5 1 28 30]); subplot(211); 
    pop_spectopo(ICA, 1, [0  ICA.pnts], 'EEG' , 'freqrange',[2 90],'electrodes','off'); 
    plot([50 50],[-100 100], 'k:');
    set(gca, 'xtick',0:10:80); title([allVhdr(n).name(8:end-5) ', Orig'])
    ICA = pop_cleanline(ICA, ...
        'bandwidth',2, ...              % 'Bandwidth',2,
        'chanlist',1:ICA.nbchan, ...    % 'ChanCompIndices',[1:ICA.nbchan],'freqrange',[1
        'computepower',0, ...           % 'ComputeSpectralPower',true,
        'linefreqs',50, ...             % 'LineFrequencies',[60 120],
        'newversion',1, ...
        'normSpectrum',0, ...           % 'NormalizeSpectrum',false,
        'p',0.01, ...                   % 'LineAlpha',0.01,
        'pad',2, ...                    % 'PaddingFactor',2,
        'plotfigures',0, ...            % 'PlotFigures',false,
        'scanforlines',1, ...           % 'ScanForLines',true,
        'sigtype','Channels', ...       % 'SignalType','Channels',
        'taperbandwidth',2, ...         % seeting for newversion
        'tau',100, ...                  % 'SmoothingFactor',100,
        'verb',0, ...                   % 'VerboseOutput',1,
        'winsize',4, ...                % 'SlidingWinLength',ICA.pnts/ICA.srate,
        'winstep',1);                   % 'SlidingWinStep',ICA.pnts/ICA.srate
    subplot(212); 
    pop_spectopo(ICA, 1, [0  ICA.pnts], 'EEG' , 'freqrange',[2 90],'electrodes','off'); 
    plot([50 50],[-100 100], 'k:');
    set(gca, 'xtick',0:10:80); title([allVhdr(n).name(8:end-5) ', cleanline newversion'])
    exportgraphics(gcf,fullfile(PATH_OUT_LN, [allVhdr(n).name(1:end-5) '_cleanline.tif']))
    close(gcf)
    
    % insert zero-channel for online-reference (according to fullRankAveRef() from Makoto Miyakoshi)
    ICA.nbchan = ICA.nbchan+1;
    ICA.data(end+1,:) = zeros(1, ICA.pnts);
    ICA.chanlocs(1,ICA.nbchan).labels = 'initialReference';
    % re-reference to average reference
    ICA = pop_reref(ICA, []);
    % discard one channel to adjust rank
    ICA = pop_select( ICA,'nochannel',{'initialReference'});
  
    % epoch data for automatic artifact rejection
    ICA = eeg_regepochs(ICA, 'recurrence', 2, 'limits', [0 2], 'rmbase', NaN);
    col=col+1; artefact_overview(1,col) = ICA.trials;  % document data set statistics
    col=col+1; artefact_overview(1,col) = ICA.pnts;    % document data set statistics
    
    % artifact rejection
    ICA = pop_eegthresh(ICA, 1, 1:ICA.nbchan ,-1000,1000,0,ICA.xmax,0,0);
    ICA = pop_jointprob(ICA, 1, 1:ICA.nbchan, 5, 2, 0, 0, 0, [], 0); % joint probability, according to Makotos recommendations and Delorme et al. 2007 (5SD; more strict for ICA)
    ICA = pop_rejkurt(  ICA, 1, 1:ICA.nbchan, 5, 2, 0, 0, 0, [], 0);          
    REJ_THRE = logical(ICA.reject.rejthresh);
    REJ_JP   = ICA.reject.rejjp;
    REJ_KURT = ICA.reject.rejkurt;
    REJ_ALL  = logical(REJ_KURT + REJ_JP + REJ_THRE);
    ICA = eeg_rejsuperpose( ICA, 1, 1, 1, 1, 1, 1, 1, 1);
    ICA = pop_rejepoch( ICA, find(ICA.reject.rejglobal), 0);    
    
    % save stats of rejection
    col=col+1; artefact_overview(1,col) = ICA.trials;
    col=col+1; artefact_overview(1,col) = artefact_overview(1,col-3) - ICA.trials;
    col=col+1; artefact_overview(1,col) = artefact_overview(1,col-2) / artefact_overview(1,col-4);
    col=col+1; artefact_overview(1,col) = 1 - artefact_overview(1,col-1);
        
    col=col+1; artefact_overview(1,col) = numel(unique([find(REJ_THRE) find(REJ_JP) find(REJ_KURT)]));
    col=col+1; artefact_overview(1,col) = sum(REJ_THRE);
    col=col+1; artefact_overview(1,col) = sum(REJ_JP);
    col=col+1; artefact_overview(1,col) = sum(REJ_KURT);
    col=col+1; artefact_overview(1,col) = numel(setdiff(find(REJ_THRE), unique([find(REJ_JP)   find(REJ_KURT)])));
    col=col+1; artefact_overview(1,col) = numel(setdiff(find(REJ_JP),   unique([find(REJ_THRE) find(REJ_KURT)])));
    col=col+1; artefact_overview(1,col) = numel(setdiff(find(REJ_KURT), unique([find(REJ_THRE) find(REJ_JP)  ])));
    col=col+1; artefact_overview(1,col) = numel(setdiff(intersect(find(REJ_THRE),find(REJ_JP)),   find(REJ_KURT)));
    col=col+1; artefact_overview(1,col) = numel(setdiff(intersect(find(REJ_THRE),find(REJ_KURT)), find(REJ_JP)));
    col=col+1; artefact_overview(1,col) = numel(setdiff(intersect(find(REJ_JP),  find(REJ_KURT)), find(REJ_THRE)));
    col=col+1; artefact_overview(1,col) = numel(intersect(intersect(find(REJ_THRE),find(REJ_JP)),find(REJ_KURT)));
    
    % set random seed for ICA
    rng('shuffle');
    seed = rng;
    col=col+1; artefact_overview(1,col) = seed.Seed;
    
    % run ICA
    tic;
    ICA = pop_runica(ICA, 'icatype', 'runica', 'extended', 1, 'interrupt', 'on', 'verbose', 'off'); % 'verbose', 'off'
    col=col+1; artefact_overview(1,col) = toc/60; % min
    
    % run ICLabel
    ICA = pop_iclabel(ICA, 'default');

    % plot IC outcome
    figure; pop_topoplot(ICA, 0, 1:ICA.nbchan, allVhdr(n).name(8:end-5),[6 6] ,0,'electrodes','on', 'iclabel', 'on');
    set(gcf, 'Units', 'centimeters', 'Position', [5 1 60 30]);
    exportgraphics(gcf, fullfile(PATH_OUT_ICA, [allVhdr(n).name(1:end-5) '_runica.tif']))
    close(gcf)
    
    % flag ICs (optional: pop_icflag(); alg. is not restricted to "winning"-ICs)
    ICA.reject.gcompreject = zeros(1, size(ICA.icaweights,1));
    classes_probs          = nan(2, size(ICA.icaweights,1));
    count_muscle = 0; count_eye = 0; count_heart = 0;
    for i=1:size(ICA.icaweights,1) % loop through all ICs
        class_prob = ICA.etc.ic_classification.ICLabel.classifications(i,:);
        max_prob   = max(class_prob);
        classes_probs(1,i) = find(class_prob == max_prob);
        classes_probs(2,i) = max_prob;
        if      (find(class_prob == max_prob) == 2 && max_prob >= .9)  || ... % find all Muscle-Winning-ICs and check if their probability is above .9
                (find(class_prob == max_prob) == 3 && max_prob >= .25) || ... % find all Eye-Winning-ICs and check if their probability is above .25
                (find(class_prob == max_prob) == 4 && max_prob >= .6)         % find all Heart-Winning-ICs and check if their probability is above .6
            ICA.reject.gcompreject(1,i) = 1;
            if find(class_prob == max_prob) == 2, count_muscle = count_muscle + 1; end
            if find(class_prob == max_prob) == 3, count_eye    = count_eye + 1;    end
            if find(class_prob == max_prob) == 4, count_heart  = count_heart + 1;  end
        end
    end
    cnt = histogram(classes_probs(1,:), 1:7+1);
    col=col+1; artefact_overview(1,col:col+6) = cnt.Values;    
    col=col+6;
    for i = 1:7 % loop over all IC classes
        col=col+1; artefact_overview(1,col) = median(classes_probs(2, classes_probs(1,:) == i));
    end
    
    % save rejection stats
    rejected_ICs = find(ICA.reject.gcompreject);
    col=col+1; artefact_overview(1,col) = count_eye;
    col=col+1; artefact_overview(1,col) = count_muscle;
    col=col+1; artefact_overview(1,col) = count_heart;
    col=col+1; artefact_overview(1,col) = numel(rejected_ICs); % document data set statistics
    
    % plot ICLabel results
    figure; pop_topoplot(ICA, 0, rejected_ICs, allVhdr(n).name(8:end-5),[] ,0,'electrodes','on', 'iclabel', 'on');
    fig_rIC = gcf;
    set(gcf, 'Position', [fig_rIC.Position(1) fig_rIC.Position(2) fig_rIC.Position(3)*1.3 fig_rIC.Position(4)])
    exportgraphics(gcf,fullfile (PATH_OUT_ICL, [allVhdr(n).name(1:end-5) '_runica_removedICs.tif']))
    close(gcf)
    
    % document data set statistics
    col=    1; artefact_epochs(1,col) = {(setdiff(find(REJ_THRE), unique([find(REJ_JP)   find(REJ_KURT)])))};
    col=col+1; artefact_epochs(1,col) = {(setdiff(find(REJ_JP),   unique([find(REJ_THRE) find(REJ_KURT)])))};
    col=col+1; artefact_epochs(1,col) = {(setdiff(find(REJ_KURT), unique([find(REJ_THRE) find(REJ_JP)  ])))};
    col=col+1; artefact_epochs(1,col) = {(setdiff(intersect(find(REJ_THRE),find(REJ_JP)),   find(REJ_KURT)))};
    col=col+1; artefact_epochs(1,col) = {(setdiff(intersect(find(REJ_THRE),find(REJ_KURT)), find(REJ_JP)))};
    col=col+1; artefact_epochs(1,col) = {(setdiff(intersect(find(REJ_JP),  find(REJ_KURT)), find(REJ_THRE)))};
    col=col+1; artefact_epochs(1,col) = {(intersect(intersect(find(REJ_THRE),find(REJ_JP)),find(REJ_KURT)))};
    col=col+1; artefact_epochs(1,col) = {unique([find(REJ_THRE) find(REJ_JP) find(REJ_KURT)])};

    % save rejection stats to hard disk
    csvwrite(fullfile(PATH_OUT_STATS, ...
        [allVhdr(n).name(1:end-5) FILE_STATS_OVERVIEW '.csv']), artefact_overview)
    out = matfile(fullfile(PATH_OUT_STATS, ...
        [allVhdr(n).name(1:end-5) FILE_STATS_EPOCHS   '.mat']), ...
        'writable',true);
    out.artefact_epochs = artefact_epochs;
    
    % ---------- Preprocessing for Time-Frequency-Analysis ----------------
    
    % filter data for final analysis
    EEG = pop_eegfiltnew(EEG, 'locutoff', 0.1); % optional: 'hicutoff'
            
    % remove noisy channels detected by first ASR
    removed_chans = setdiff({EEG.chanlocs.labels}, {ICA.chanlocs.labels});
    EEG = pop_select(EEG, 'rmchannel', removed_chans);
    
    % remove line noise
    EEG = pop_cleanline(EEG, ...
        'bandwidth',2, ...              % 'Bandwidth',2,
        'chanlist',1:EEG.nbchan, ...    % 'ChanCompIndices',[1:EEG.nbchan],'freqrange',[1
        'computepower',0, ...           % 'ComputeSpectralPower',true,
        'linefreqs',50, ...             % 'LineFrequencies',[60 120],
        'newversion',1, ...
        'normSpectrum',0, ...           % 'NormalizeSpectrum',false,
        'p',0.01, ...                   % 'LineAlpha',0.01,
        'pad',2, ...                    % 'PaddingFactor',2,
        'plotfigures',0, ...            % 'PlotFigures',false,
        'scanforlines',1, ...           % 'ScanForLines',true,
        'sigtype','Channels', ...       % 'SignalType','Channels',
        'taperbandwidth',2, ...         % setting for newversion
        'tau',100, ...                  % 'SmoothingFactor',100,
        'verb',0, ...                   % 'VerboseOutput',1,
        'winsize',4, ...                % 'SlidingWinLength',EEG.pnts/EEG.srate,
        'winstep',1);                   % 'SlidingWinStep',EEG.pnts/EEG.srate
    
    % insert zero-channel for online-reference (according to fullRankAveRef() from Makoto Miyakoshi)
    EEG.nbchan = EEG.nbchan+1;
    EEG.data(end+1,:) = zeros(1, EEG.pnts);
    EEG.chanlocs(1,EEG.nbchan).labels = 'initialReference';
    % re-reference to average reference
    EEG = pop_reref(EEG, []);
    % discard one channel to adjust rank
    EEG = pop_select( EEG,'nochannel',{'initialReference'});
    
    % copy IC weights to EEG data set
    EEG.icachansind              = ICA.icachansind;
    EEG.icaweights               = ICA.icaweights;
    EEG.icasphere                = ICA.icasphere;
    EEG.icawinv                  = ICA.icawinv;
    EEG.etc.ic_classification    = ICA.etc.ic_classification;    
    EEG.etc.icaweights_beforerms = ICA.etc.icaweights_beforerms;
    EEG.etc.icasphere_beforerms  = ICA.etc.icasphere_beforerms;
    EEG.reject.gcompreject       = ICA.reject.gcompreject;
    
    % remove components previously marked for rejection
    EEG = pop_subcomp( EEG, [], 0); % remove IC saved in EEG.reject.gcompreject
    
    % interpolate channels (removed and online reference)
    EEG = pop_interp(EEG, orig_chanlocs, 'spherical');
    
    % save continuous preprocessed EEG dataset
    saveName = [allVhdr(n).name(1:end-5) FILE_OUT];
    EEG.setname = saveName;
    EEG = pop_saveset(EEG, 'filename', saveName, 'filepath', PATH_OUT_PREPRO);
   
end

