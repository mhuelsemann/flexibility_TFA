%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Title:         mass-univariate analysis with cluster-based permutation  %
%                tests of time-frequency data                             %
% Subtitle:      save empirical data for reliability                      %
% Author:        Mareike J. Hülsemann, University of Mainz                %
% Date:          01/2024                                                  %
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

% define PATHs
[filepath,~,~] = fileparts(matlab.desktop.editor.getActiveFilename);
parts          = regexp(filepath,filesep,'split');
PATH_MAIN      = fullfile(parts{1:end-1});
PATH_MAIN      = fullfile(PATH_MAIN, 'derivatives', task_name, ref);

PATH_IN   = fullfile(PATH_MAIN, 'tf_power');
PATH_OUT  = fullfile(PATH_MAIN, 'tf_stats');


%% ---------- Load Data ---------------------------------------------------

% list all data sets
allData = dir(fullfile(PATH_IN, '*.mat'));
% exclude data sets
allData = allData(~contains({allData.name}, subj2excl)); 

% load data
load(fullfile(PATH_IN, allData(1).name), 'chanlocs', 'frex', 'trial_time','trial_info')
all_eegpower = nan(numel(allData),numel(chanlocs),numel(frex),numel(trial_time),numel(trial_info));
for n=1:numel(allData)
    load(fullfile(PATH_IN, allData(n).name), 'eegpower');
    all_eegpower(n,:,:,:,:) = eegpower;
end
clear eegpower


%% ---------- Save Output -------------------------------------------------

save(fullfile(PATH_OUT, [allData(1).name(1:6) task_name '_tf_power.mat']), ...
    'all_eegpower', 'allData', 'note', 'task_name', 'ref', ...
    'chanlocs', 'frex', 'trial_time', 'trial_info');


