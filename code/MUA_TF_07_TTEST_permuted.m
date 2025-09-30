%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Title:         mass-univariate analysis with cluster-based permutation  %
%                tests of time-frequency data                             %
% Subtitle:      t-tests for permuted data                                %
% Author:        Mareike J. Hülsemann, University of Mainz                %
% Date:          12/2023                                                  %
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

% definitions
task_name     = 'oddeven'; % 'oddeven', 'globallocal', 'numberletter'
ref           = 'cue'; % 'cue', 'target'
switch task_name
    case 'oddeven'
        subj2excl = {'0022'; '0090'; '0121'};
    case 'globallocal'
        subj2excl = {'0013'; '0090'};
    case 'numberletter'
        subj2excl = {'0010'; '0054'; '0106'};
    otherwise
        warning('Typo in variable [task_name].')
end
note      = 'same amount of trials';
nperm     = 100; % 50, 100
nr        = 9; % 0-9

% define PATHs
[filepath,~,~] = fileparts(matlab.desktop.editor.getActiveFilename);
parts          = regexp(filepath,filesep,'split');
PATH_MAIN      = fullfile(parts{1:end-1});
PATH_MAIN      = fullfile(PATH_MAIN, 'derivatives', task_name, ref);

PATH_IN   = fullfile(PATH_MAIN, 'tf_power_perm');
PATH_OUT  = fullfile(PATH_MAIN, 'tf_stats');

% define FILE names
FILE_OUT = ['Exp23_' task_name '_tf_power_TTEST_perm_' num2str(nr) '.mat'];


%% ---------- Load Data ---------------------------------------------------

% list all data sets
allData = dir(fullfile(PATH_IN, ['*' num2str(nr) '.mat']));
% exclude data sets
allData = allData(~contains({allData.name}, subj2excl));

% load data
load(fullfile(PATH_IN, allData(1).name), 'chanlocs', 'frex', 'trial_time')
all_eegpower = nan(numel(allData),numel(chanlocs),numel(frex),numel(trial_time),2,nperm);
for n=1:numel(allData)
    load(fullfile(PATH_IN, allData(n).name), 'eegpower_perm');
    all_eegpower(n,:,:,:,:,:) = eegpower_perm;
    disp(['Loaded ' num2str(n) ' from ' num2str(numel(allData)) '.'])
end
clear eegpower_perm


%% ---------- Execute t-Tests ---------------------------------------------

% initialise
TTEST_perm_vars = {'T', 'df', 'p'};
TTEST_perm      = nan(numel(chanlocs), numel(frex), numel(trial_time), numel(TTEST_perm_vars), nperm);

for pi=1:nperm
    parfor ci=1:numel(chanlocs)
        tmp_ttest = nan(numel(frex), numel(trial_time), numel(TTEST_perm_vars));
        data = squeeze(all_eegpower(:,ci,:,:,:,pi));

        for fi=1:numel(frex)
            for ti=1:numel(trial_time)

                % get data
                data1 = squeeze(data(:,fi,ti,1));
                data2 = squeeze(data(:,fi,ti,2));

                % run t-tests
                [H,P,CI,STATS] = ttest(data1, data2);

                % save stats
                col =     1; tmp_ttest(fi,ti,col) = STATS.tstat;
                col = col+1; tmp_ttest(fi,ti,col) = STATS.df;
                col = col+1; tmp_ttest(fi,ti,col) = P;

            end
        end
        TTEST_perm(ci,:,:,:,pi) = tmp_ttest;
    end
    disp(['Permutation ' num2str(pi) '; ' datestr(datetime('now')) '.'])
end


%% ---------- save output -------------------------------------------------

% infos
note2 = 'channels x frequencies x time x test statistics x nperm';

% save output
save(fullfile(PATH_OUT, FILE_OUT), 'allData', 'TTEST_perm', 'TTEST_perm_vars', ...
    'note', 'note2', 'chanlocs', 'frex', 'trial_time');



%% ---------- combine output ----------------------------------------------
tmp = nan(numel(chanlocs), numel(frex), numel(trial_time), numel(TTEST_perm_vars), 1000);

for i=0:9
    load(fullfile(PATH_OUT, [FILE_OUT(1:end-5) num2str(i) '.mat']))
    tmp(:,:,:,:,i*nperm+1:(i+1)*nperm) = TTEST_perm;
end
TTEST_perm = tmp;

% save output
save(fullfile(PATH_OUT, [FILE_OUT(1:end-6) '.mat']), 'allData', ...
    'TTEST_perm', 'TTEST_perm_vars', ...
    'note', 'note2', 'chanlocs', 'frex', 'trial_time', ...
    '-v7.3');
