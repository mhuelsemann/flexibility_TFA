%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Title:         mass-univariate analysis with cluster-based permutation  %
%                tests of time-frequency data                             %
% Subtitle:      wavelets                                                 %
% Author:        Mareike J. Hülsemann, University of Mainz                %
% Date:          11/2023                                                  %
% Last Modified: 09/2025                                                  %
% Note:          Exp23, Switching Tasks                                   %
%                Original Wavelet Code by Mike X Cohen (code accompanying %
%                the book "Analyzing Neural Time Series Data", chapter 13)%
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
ref           = 'cue'; % 'cue', 'target'
save_rawpower = 1; % choose with 0/1 whether to save large rawpower file
switch task_name
    case 'oddeven'
        marker_start  = 'S  6';
        letters2remove = 36;
    case 'globallocal'        
        marker_start  = 'S 14';
        letters2remove = 40;
    case 'numberletter'
        marker_start  = 'S 12';
        letters2remove = 43;
    otherwise
        warning('Typo in variable [task_name].')
end
marker_start_num = str2double(marker_start(3:end));
switch ref
    case 'cue'
        event_indices = {'S 31' 'S 32' 'S 33' 'S 34'};
        epoch_limits    = [-1.75 2.0]; % sec
        trial_limits    = [-750 1000]; % ms
        baseline_limits = [-750 -250]; % ms
    case 'target'
        event_indices = {'S 51' 'S 52' 'S 53' 'S 54'};
        epoch_limits    = [-2.95 2.5]; % sec
        trial_limits    = [-1950 1500]; % ms
        baseline_limits = [-1950 -1450]; % ms
    otherwise
        warning('Typo in variable [ref].')
end
trial_info = {'switch,all', ...
        'repeat,all', ...
        'switch,odd', ...
        'repeat,odd', ...
        'switch,even', ...
        'repeat,even'};
% IMPORTANT NOTE: In this script and in output files, condition labels of
% 'switch_even' and 'repeat_odd' have been mixed up (see line 249-250).

% define PATHs
[filepath,~,~] = fileparts(matlab.desktop.editor.getActiveFilename);
parts          = regexp(filepath,filesep,'split');
PATH_MAIN      = fullfile(parts{1:end-1});

PATH_EEGLAB    = [PATH_MAIN '/eeglab2023.0/'];
PATH_HD        = PATH_MAIN;

PATH_MAIN      = fullfile(PATH_MAIN, 'derivatives', task_name, ref);
PATH_IN        = fullfile(PATH_MAIN, 'preprocessed');
PATH_OUT       = fullfile(PATH_MAIN, 'tf_power');
PATH_OUT_STATS = fullfile(PATH_MAIN, 'stats');
PATH_OUT_RAW   = fullfile(PATH_HD, 'derivatives', task_name, ref, 'tf_rawpower');

% define FILE names
FILE_MAIN = ['_epoched_' ref '_' num2str(trial_limits(1)) '_' num2str(trial_limits(2))];


% load list of preprocessed EEG data
allSets = dir(fullfile(PATH_IN, '*.set'));

% initialize preprocessing stats
artefact_varnames = {...
    'id', ...
    'task', ...
    'srate_for_TF', ...
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
    'nb_trials_rejected_inTPK'};
artefact_epochs_varnames = {'trials_rejected_unique_thre', ...
    'trials_rejected_unique_prob', ...
    'trials_rejected_unique_kurt', ...
    'trials_rejected_inTP_notK', ...
    'trials_rejected_inTK_notP', ...
    'trials_rejected_inPK_notT', ...
    'trials_rejected_inTPK', ...
    'trials_rejected_all'};

% wavelet parameter
min_freq   =  3;
max_freq   = 60;
num_frex   = 25;
min_cyc    =  3;
max_cyc    =  6;


%% ---------- start EEGLAB ------------------------------------------------
addpath(PATH_EEGLAB);
eeglab; close(gcf);
pop_editoptions('option_single', 0, ...
    'option_savetwofiles', 0, ... % new default: FALSE
    'option_computeica',   1, ... % one only has EEG.icaact if TRUE
    'option_checkversion', 0);


%% ---------- Loop Over Subjects ------------------------------------------
for n=1:length(allSets)
    disp(['**************** Started Set ' num2str(n) ' of ' num2str(length(allSets)) '; ' allSets(n).name '; ' datestr(datetime()) ' ****************'])
    
    % initialise
    artefact_overview  = nan(1,  numel(artefact_varnames));
    artefact_epochs    = cell(1, numel(artefact_epochs_varnames));

    % initialise
    frex     = logspace(log10(min_freq),log10(max_freq),num_frex);
    s        = logspace(log10(min_cyc), log10(max_cyc), num_frex)./(2*pi*frex); % SD of Gaussian
    
    % ---------- Loading data  --------------------------------------------
    % load data set
    EEG = pop_loadset('filename',allSets(n).name,'filepath',PATH_IN);

    
    % ---------- Preprocessing for Time-Frequency-Analysis ----------------
    
    % save stats of rejection
    col=    1; artefact_overview(1,col) = str2double(allSets(n).name(7:end-letters2remove));
    col=col+1; artefact_overview(1,col) = marker_start_num;
    col=col+1; artefact_overview(1,col) = EEG.srate;
    
    % extract trials (no baseline removal)
    EEG = pop_epoch(EEG, event_indices, epoch_limits);
    if strcmp(task_name,'oddeven') && n==2
        % remove last trial for subject 0002, as last trial is incomplete (no RT)
        EEG = pop_select(EEG, 'trial', 1:267);
    end 
    
    % detect artefactual epochs
    EEG = pop_eegthresh(EEG, 1, 1:EEG.nbchan ,-1000,1000,0,EEG.xmax,0,0);
    EEG = pop_jointprob(EEG, 1, 1:EEG.nbchan, 5, 2, 0, 0, 0, [], 0); % joint probability, according to Makotos recommendations and Delorme et al. 2007 (5SD; more strict for EEG)
    EEG = pop_rejkurt(  EEG, 1, 1:EEG.nbchan, 5, 2, 0, 0, 0, [], 0);          
    REJ_THRE = logical(EEG.reject.rejthresh);
    REJ_JP   = EEG.reject.rejjp;
    REJ_KURT = EEG.reject.rejkurt;
    REJ_ALL  = logical(REJ_KURT + REJ_JP + REJ_THRE);
    EEG = eeg_rejsuperpose( EEG, 1, 1, 1, 1, 1, 1, 1, 1);
        
    % save artefactual epochs
    if strcmp(task_name, 'numberletter') && strcmp(ref, 'target') && ismember(n, [56, 66, 72])
        EEG.etc.epoch_info(end,:) = [];
    elseif strcmp(task_name, 'numberletter') && strcmp(ref, 'target') && ismember(n, 110)
        EEG.etc.epoch_info(1,:) = [];
    end
    EEG.etc.epoch_info(:,10) = EEG.reject.rejglobal;
    
    % save stats of rejection
    col=col+1; artefact_overview(1,col) = EEG.trials;
    col=col+1; artefact_overview(1,col) = EEG.pnts;
    col=col+1; artefact_overview(1,col) = EEG.trials - sum(EEG.reject.rejglobal);
    col=col+1; artefact_overview(1,col) = sum(EEG.reject.rejglobal);
    col=col+1; artefact_overview(1,col) = 1 - sum(EEG.reject.rejglobal)/EEG.trials;
    col=col+1; artefact_overview(1,col) = sum(EEG.reject.rejglobal)/EEG.trials;
        
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
       
    % document data set statistics
    col=    1; artefact_epochs(1,col) = {(setdiff(find(REJ_THRE), unique([find(REJ_JP)   find(REJ_KURT)])))};
    col=col+1; artefact_epochs(1,col) = {(setdiff(find(REJ_JP),   unique([find(REJ_THRE) find(REJ_KURT)])))};
    col=col+1; artefact_epochs(1,col) = {(setdiff(find(REJ_KURT), unique([find(REJ_THRE) find(REJ_JP)  ])))};
    col=col+1; artefact_epochs(1,col) = {(setdiff(intersect(find(REJ_THRE),find(REJ_JP)),   find(REJ_KURT)))};
    col=col+1; artefact_epochs(1,col) = {(setdiff(intersect(find(REJ_THRE),find(REJ_KURT)), find(REJ_JP)))};
    col=col+1; artefact_epochs(1,col) = {(setdiff(intersect(find(REJ_JP),  find(REJ_KURT)), find(REJ_THRE)))};
    col=col+1; artefact_epochs(1,col) = {(intersect(intersect(find(REJ_THRE),find(REJ_JP)),find(REJ_KURT)))};
    col=col+1; artefact_epochs(1,col) = {unique([find(REJ_THRE) find(REJ_JP) find(REJ_KURT)])};    

    % save preprocessing stats
    csvwrite(fullfile(PATH_OUT_STATS, [allSets(n).name(1:end-letters2remove) FILE_MAIN '_artefact_overview.csv']), artefact_overview)
    out1 = matfile(fullfile(PATH_OUT_STATS, [allSets(n).name(1:end-letters2remove) FILE_MAIN '_artefact_epochs.mat']), 'writable', true);
    out1.artefact_epochs          = artefact_epochs;
    out1.artefact_epochs_varnames = artefact_epochs_varnames;
    
    % reject trials with EEG artifacts
    EEG.etc.epoch_info(EEG.reject.rejglobal,:) = [];
    EEG = pop_rejepoch(EEG, find(EEG.reject.rejglobal), 0);
    
    % reject trials with incorrect responses
    EEG = pop_rejepoch(EEG, find(EEG.etc.epoch_info(:, 5) == 0), 0);
    EEG.etc.epoch_info(EEG.etc.epoch_info(:, 5) == 0,:) = [];
    
    % reject trials with RT outliers
    EEG = pop_rejepoch(EEG, find(EEG.etc.epoch_info(:, 8) == 1), 0);
    EEG.etc.epoch_info(EEG.etc.epoch_info(:, 8) == 1,:) = [];
    
    
    % ---------- Define trials for conditions ---------------------
    
    tmp   = cell(4,2);
    trial = cell(numel(trial_info),1);

    tmp{1,1} = find(EEG.etc.epoch_info(:,2) == 1 & EEG.etc.epoch_info(:,3) == 1); % Switch, LessMore
    tmp{2,1} = find(EEG.etc.epoch_info(:,2) == 1 & EEG.etc.epoch_info(:,3) == 0); % Repeat, LessMore
    tmp{3,1} = find(EEG.etc.epoch_info(:,2) == 2 & EEG.etc.epoch_info(:,3) == 1); % Switch, OddEven
    tmp{4,1} = find(EEG.etc.epoch_info(:,2) == 2 & EEG.etc.epoch_info(:,3) == 0); % Repeat, OddEven

    for i=1:size(tmp,1)
        tmp{i,2} = numel(tmp{i,1});
    end
    lb = min([tmp{:,2}]);

    for i=1:size(tmp,1)
        selection = tmp{i,1};
        selection = selection(randperm(numel(tmp{i,1})));
        tmp{i,1}  = selection(1:lb);
    end
    
    trial{1,1} = sort([tmp{1,1}; tmp{3,1}]); % Switch
    trial{2,1} = sort([tmp{2,1}; tmp{4,1}]); % Repeat
    trial{3,1} = trial{1,1}(1:2:end);  % Switch, reliability
    trial{4,1} = trial{1,1}(2:2:end);  % Switch, reliability
    trial{5,1} = trial{2,1}(1:2:end);  % Repeat, reliability
    trial{6,1} = trial{2,1}(2:2:end);  % Repeat, reliability
    
        
    % ---------- Time-Frequency-Analysis via Wavelets ---------------------
    
    % settings
    pnts       = EEG.pnts;
    trials     = EEG.trials;
    trialidx1  = dsearchn(EEG.times',trial_limits(1));
    trialidx2  = dsearchn(EEG.times',trial_limits(2));
    trial_time = EEG.times(trialidx1:trialidx2);
    baseidx1   = dsearchn(trial_time',baseline_limits(1));
    baseidx2   = dsearchn(trial_time',baseline_limits(2))-1;
    
    % define wavelet parameters
    rd_pnts = 1000/EEG.srate; % fill up timeepoch'
    time_1  = EEG.times(1)/EEG.srate;
    time_2  = (EEG.times(end)+rd_pnts)/EEG.srate;
    time    = time_1:1/EEG.srate:time_2;
            
    % definte convolution parameters
    n_wavelet            = length(time);
    n_data               = EEG.pnts*EEG.trials;
    n_convolution        = n_wavelet+n_data-1;
    n_conv_pow2          = pow2(nextpow2(n_convolution));
    half_of_wavelet_size = (n_wavelet-1)/2;
    
    % initialize
    rawpower    = nan(EEG.nbchan, num_frex, trialidx2-trialidx1+1, EEG.trials);        % channels X frequencies X time X trials
    eegpower    = nan(EEG.nbchan, num_frex, trialidx2-trialidx1+1, numel(trial_info)); % channels X frequencies X time X [conditions * allTrials/odd/even]
    
    % loop through all channels
    for ci=1:EEG.nbchan
        % get FFT of data
        eegfft = fft(reshape(EEG.data(ci,:,:),1,EEG.pnts*EEG.trials),n_conv_pow2);
        
        % loop through frequencies
        parfor fi=1:num_frex
            % build wavelet
            wavelet = fft( exp(2*1i*pi*frex(fi).*time) .* exp(-time.^2./(2*(s(fi)^2))) , n_conv_pow2 );  % formula 13.10, p. 157 without scaling factor (not needed when applying baseline correction)
            
            % convolution
            eegconv = ifft(wavelet.*eegfft);
            eegconv = eegconv(1:n_convolution);
            eegconv = eegconv(half_of_wavelet_size+1:end-half_of_wavelet_size);
            eegconv = reshape(eegconv,pnts,trials);
            
            % remove edges of time window (=buffer time)
            rawpower(ci,fi,:,:) = eegconv(trialidx1:trialidx2,:);            
        end
    end
    
    % extract power
    rawpower = abs(rawpower).^2;
    
    % save power for later permutations
    if save_rawpower
        out2 = matfile(fullfile(PATH_OUT_RAW, [allSets(n).name(1:end-letters2remove) FILE_MAIN '_tf_rawpower.mat']), 'writable', true);
        out2.rawpower   = rawpower;
        out2.task       = task_name;
        out2.info       = ['channels X frequencies X time X epochs; ' ref ' = 0'];
        out2.chanlocs   = EEG.chanlocs;
        out2.frex       = frex;
        out2.s          = s;
        out2.trial_time = trial_time;
        out2.epoch_info = EEG.etc.epoch_info;
        out2.epoch_info_varnames = EEG.etc.epoch_info_varnames;
        out2.trial_info = trial_info;
        out2.trial      = trial;
    end
    
    % calculate mean for baseline over all included trials (condition average)
    baseline = mean( mean( rawpower(:, :, baseidx1:baseidx2, sort([trial{1,1}; trial{2,1}]) ) ,4) ,3);
    % 1. mean over 4th dimension = mean over all trials
    % 2. mean over 3rd dimension = mean over the specified baseline period
    % => results in a 2-dim matrix: channels x frequencies
    
    % loop over all trial types
    for i = 1:size(trial,1)
        % average power over trials (no baseline transform)
        eegpower(:,:,:,i) = mean( rawpower(:,:,:,trial{i,1}) ,4);
        
        % baseline transform - condition average
        eegpower(:,:,:,i)  = 10*log10( eegpower(:,:,:,i) ./ baseline );
        % 3-dim matrix is point-wise divided by 2-dim matrix; 
        % first two dimension correspond; 2-dim matrix is replicated in the
        % 3rd dimension.        
    end
    
    % downsample data
    ds_rate              = 20; % 10ms ~ 100 Hz, 20ms ~ 50 Hz, 25ms ~ 40 Hz (lowest)
    ds_timesidx          = dsearchn(trial_time',(trial_limits(1):ds_rate:trial_limits(2))');
    eegpower             = eegpower(:,:,ds_timesidx,:);
    trial_time           = trial_time(ds_timesidx);
    
    % save TF data
    out3 = matfile(fullfile(PATH_OUT, [allSets(n).name(1:end-letters2remove) FILE_MAIN '_tf_power.mat']), 'writable', true);
    out3.eegpower             = eegpower;
    out3.task                 = task_name;
    out3.info                 = ['channels X frequencies X time X trialtypes; ' ref ' = 0'];
    out3.chanlocs             = EEG.chanlocs;
    out3.frex                 = frex;
    out3.s                    = s;
    out3.trial_time           = trial_time;
    out3.trial_info           = trial_info;
    out3.trial                = trial;
end

