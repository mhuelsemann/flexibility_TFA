%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Title:         mass-univariate analysis with cluster-based permutation  %
%                tests of time-frequency data                             %
% Subtitle:      wavelets on permuted data                                %
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
run           = 0;
nperm         = 100;

% define PATHs
[filepath,~,~] = fileparts(matlab.desktop.editor.getActiveFilename);
parts          = regexp(filepath,filesep,'split');
PATH_MAIN      = fullfile(parts{1:end-1});
PATH_MAIN      = fullfile(PATH_MAIN, 'derivatives', task_name, ref);

PATH_IN        = fullfile(PATH_MAIN, 'tf_rawpower');
PATH_OUT       = fullfile(PATH_MAIN, 'tf_power_perm');

% define FILE names
FILE_OUT       = ['_tfp_perm_' num2str(run)];


% load list of preprocessed EEG data
allMats = dir(fullfile(PATH_IN, '*.mat'));

% downsampling data
switch ref
    case 'cue'
        trial_limits    = [-750 1000]; % ms
        baseline_limits = [-750 -250]; % ms
    case 'target'
        trial_limits    = [-1950  1500]; % ms
        baseline_limits = [-1950 -1450]; % ms
    otherwise
        warning('Typo in variable [ref].')
end
sr              = 1000;
trial_time      = trial_limits(1):1000/sr:trial_limits(2); % 1000 for ms
ds_rate         = 20; % 10ms ~ 100 Hz, 20ms ~ 50 Hz, 25ms ~ 40 Hz (lowest)
ds_timesidx     = dsearchn(trial_time',(trial_limits(1):ds_rate:trial_limits(2))');
baseidx1        = dsearchn(trial_time',baseline_limits(1));
baseidx2        = dsearchn(trial_time',baseline_limits(2))-1;


%% ---------- Loop Over Subjects ------------------------------------------

for n=1:length(allMats)
% parfor n=1:length(allMats)
    
    % ---------- shuffle RNG ----------------------------------------------
    rng('shuffle');
    g    = rng;
    seed = g.Seed;
   
    
    % ---------- Loading data  --------------------------------------------
    load(fullfile(PATH_IN, allMats(n).name))
    %[rawpower, trial, epoch_info, trial_info, trial_time, chanlocs, frex, s] = parload_rawpower(fullfile(PATH_IN, allMats(n).name));
    
    
    % ---------- Permuting data -------------------------------------------
    
    % initialise
    eegpower_perm             = nan(size(rawpower,1), size(rawpower,2), numel(ds_timesidx), 2, nperm); % channels X frequencies X time X conditions X permutations
    
    % get all trials used in original data
    trials1 = trial{1,1};
    trials2 = trial{2,1};    
    ntrial = numel([trials1; trials2]);
    
%     parfor pi = 1:nperm
    for pi = 1:nperm
    	
        % initialize
        trials                   = sort([trials1; trials2]);
        tmp_eegpower             = nan(size(rawpower,1), size(rawpower,2), size(rawpower,3), 2); % channels X frequencies X time X conditions
        
        % randomize trial order
        trials = trials(randperm(ntrial));
        
        % average power over trials (no baseline transform)
        tmp_eegpower(:,:,:,1) = mean(rawpower(:,:,:,trials(1:ntrial/2)    ),4);
        tmp_eegpower(:,:,:,2) = mean(rawpower(:,:,:,trials(ntrial/2+1:end)),4);
        
        % calculate mean for baseline over all trials (condition average)
        baseline = mean( mean( rawpower(:,:,baseidx1:baseidx2,trials) ,4) ,3);
        % 1. mean over 4th dimension = mean over all trials
        % 2. mean over 3rd dimension = mean over the specified baseline period
        % => results in a 2-dim matrix: channels x frequencies
            
        % baseline transform - condition average
        tmp_eegpower(:,:,:,1)  = 10*log10( tmp_eegpower(:,:,:,1) ./ baseline );
        tmp_eegpower(:,:,:,2)  = 10*log10( tmp_eegpower(:,:,:,2) ./ baseline );
        % 3-dim matrix is point-wise divided by 2-dim matrix;
        % first two dimension correspond; 2-dim matrix is replicated in the
        % 3rd dimension.
            
        
        % downsample data
        eegpower_perm(:,:,:,:,pi)             = tmp_eegpower(:,:,ds_timesidx,:);
       
%         disp(pi)
    end

    
    % save permuted TF data
    out = matfile(fullfile(PATH_OUT, [allMats(n).name(1:end-16) FILE_OUT '.mat']), 'writable', true);
    out.eegpower_perm             = eegpower_perm;
    out.task                      = task_name;
    out.info                      = ['channels X frequencies X time X trialtypes X permutations; ' ref ' = 0'];
    out.chanlocs                  = chanlocs;
    out.frex                      = frex;
    out.s                         = s;
    out.trial_time                = trial_time(ds_timesidx);
    out.trial_info                = trial_info(1:2);
    out.trial                     = trial;
    out.seed                      = seed;
    
    disp(['################ Finished Set ' num2str(n) ': ' allMats(n).name ' ################'])
    
end



%% ---------- internal functions ------------------------------------------
% function [rawpower, trial, epoch_info, trial_info, trial_time, ...
%     chanlocs, frex, s] = parload_rawpower(file)
% S = load( file );
% rawpower = S.rawpower;
% trial = S.trial;
% epoch_info = S.epoch_info;
% trial_info = S.trial_info;
% trial_time = S.trial_time;
% chanlocs   = S.chanlocs;
% frex       = S.frex;
% s          = S.s;
% end



