%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Title:         mass-univariate analysis with cluster-based permutation  %
%                tests of time-frequency data                             %
% Subtitle:      export data to R                                         %
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
task_name = 'numberletter'; % 'oddeven', 'globallocal', 'numberletter'
ref = 'target'; % 'cue', 'target'
cl_name   = 'Beta'; % 'Theta', 'Alpha', 'Beta'

p = '.005';     % for initial univariate t-tests
OPTION = '1SD'; % for core region

% define PATHs
[filepath,~,~] = fileparts(matlab.desktop.editor.getActiveFilename);
parts          = regexp(filepath,filesep,'split');
PATH_MAIN      = fullfile(parts{1:end-1});

PATH_IN   = fullfile(PATH_MAIN, 'derivatives', task_name, ref, 'tf_stats');

% define FILEs
FILE_IN_POWR  = ['Exp23_' task_name '_tf_power.mat'];
FILE_IN_CLST  = ['Exp23_' task_name '_tf_power_cluster_p' p(2:end) '.mat'];
FILE_OUT_POWR = ['Exp23_' task_name '_tf_power_core' OPTION '_' cl_name '.csv'];


%% ---------- Summarise TF Power ------------------------------------------

% load power data
load(fullfile(PATH_IN, FILE_IN_POWR), 'all_eegpower', 'allData', ...
    'chanlocs','frex','trial_time')

% load cluster data
load(fullfile(PATH_IN, FILE_IN_CLST), 'cluster')


% define subclusters to use
if     strcmp(task_name, 'oddeven')      && strcmp(ref, 'cue')    && strcmp(cl_name, 'Theta')
    subclst = 1;
    % sanity check
    if cluster.subcluster(subclst).clsize(cluster.subcluster(subclst).largest) ~= 27 || ...
            cluster.subcluster(subclst).tmass(cluster.subcluster(subclst).largest) < 98.90 || ...
            cluster.subcluster(subclst).tmass(cluster.subcluster(subclst).largest) > 98.92
        warning('Check whether correct cluster was selected.')
    end
elseif strcmp(task_name, 'globallocal')  && strcmp(ref, 'cue')    && strcmp(cl_name, 'Theta')
    subclst = 1;
    % sanity check
    if cluster.subcluster(subclst).clsize(cluster.subcluster(subclst).largest) ~= 56 || ...
            cluster.subcluster(subclst).tmass(cluster.subcluster(subclst).largest) < 229.47 || ...
            cluster.subcluster(subclst).tmass(cluster.subcluster(subclst).largest) > 229.49
        warning('Check whether correct cluster was selected.')
    end
elseif strcmp(task_name, 'numberletter') && strcmp(ref, 'cue')    && strcmp(cl_name, 'Theta')
    subclst = 1;
    % sanity check
    if cluster.subcluster(subclst).clsize(cluster.subcluster(subclst).largest) ~= 21 || ...
            cluster.subcluster(subclst).tmass(cluster.subcluster(subclst).largest) < 82.02 || ...
            cluster.subcluster(subclst).tmass(cluster.subcluster(subclst).largest) > 82.04
        warning('Check whether correct cluster was selected.')
    end

elseif strcmp(task_name, 'oddeven')      && strcmp(ref, 'cue')    && strcmp(cl_name, 'Alpha')
    subclst = 20;
    % sanity check
    if cluster.subcluster(subclst).clsize(cluster.subcluster(subclst).largest) ~= 1047 || ...
            cluster.subcluster(subclst).tmass(cluster.subcluster(subclst).largest) > -6590.42 || ...
            cluster.subcluster(subclst).tmass(cluster.subcluster(subclst).largest) < -6590.44
        warning('Check whether correct cluster was selected.')
    end
elseif strcmp(task_name, 'globallocal')  && strcmp(ref, 'cue')    && strcmp(cl_name, 'Alpha')
    subclst = 29;
    % sanity check
    if cluster.subcluster(subclst).clsize(cluster.subcluster(subclst).largest) ~= 1640 || ...
            cluster.subcluster(subclst).tmass(cluster.subcluster(subclst).largest) > -11250.79 || ...
            cluster.subcluster(subclst).tmass(cluster.subcluster(subclst).largest) < -11250.81
        warning('Check whether correct cluster was selected.')
    end
elseif strcmp(task_name, 'numberletter') && strcmp(ref, 'cue')    && strcmp(cl_name, 'Alpha')
    subclst = 15;
    % sanity check
    if cluster.subcluster(subclst).clsize(cluster.subcluster(subclst).largest) ~= 1028 || ...
            cluster.subcluster(subclst).tmass(cluster.subcluster(subclst).largest) > -6252.65 || ...
            cluster.subcluster(subclst).tmass(cluster.subcluster(subclst).largest) < -6252.67
        warning('Check whether correct cluster was selected.')
    end

elseif strcmp(task_name, 'oddeven')      && strcmp(ref, 'target') && strcmp(cl_name, 'Theta')
    subclst = 12;
    % sanity check
    if cluster.subcluster(subclst).clsize(cluster.subcluster(subclst).largest) ~= 62 || ...
            cluster.subcluster(subclst).tmass(cluster.subcluster(subclst).largest) < 283.05 || ...
            cluster.subcluster(subclst).tmass(cluster.subcluster(subclst).largest) > 283.06
        warning('Check whether correct cluster was selected.')
    end
elseif strcmp(task_name, 'numberletter') && strcmp(ref, 'target') && strcmp(cl_name, 'Theta')
    subclst = 1;
    % sanity check
    if cluster.subcluster(subclst).clsize(cluster.subcluster(subclst).largest) ~= 51 || ...
            cluster.subcluster(subclst).tmass(cluster.subcluster(subclst).largest) < 200.14 || ...
            cluster.subcluster(subclst).tmass(cluster.subcluster(subclst).largest) > 200.15
        warning('Check whether correct cluster was selected.')
    end
elseif strcmp(task_name, 'oddeven')      && strcmp(ref, 'target') && strcmp(cl_name, 'Alpha')
    subclst = 49;
    % sanity check
    if cluster.subcluster(subclst).clsize(cluster.subcluster(subclst).largest) ~= 1293 || ...
            cluster.subcluster(subclst).tmass(cluster.subcluster(subclst).largest) > -7020.32 || ...
            cluster.subcluster(subclst).tmass(cluster.subcluster(subclst).largest) < -7020.33
        warning('Check whether correct cluster was selected.')
    end
elseif strcmp(task_name, 'globallocal')  && strcmp(ref, 'target') && strcmp(cl_name, 'Alpha')
    subclst = 46;
    % sanity check
    if cluster.subcluster(subclst).clsize(cluster.subcluster(subclst).largest) ~= 2432 || ...
            cluster.subcluster(subclst).tmass(cluster.subcluster(subclst).largest) > -15431.96 || ...
            cluster.subcluster(subclst).tmass(cluster.subcluster(subclst).largest) < -15431.98
        warning('Check whether correct cluster was selected.')
    end
elseif strcmp(task_name, 'numberletter') && strcmp(ref, 'target') && strcmp(cl_name, 'Alpha')
    subclst = 21;
    % sanity check
    if cluster.subcluster(subclst).clsize(cluster.subcluster(subclst).largest) ~= 1641 || ...
            cluster.subcluster(subclst).tmass(cluster.subcluster(subclst).largest) > -8914.79 || ...
            cluster.subcluster(subclst).tmass(cluster.subcluster(subclst).largest) < -8914.81
        warning('Check whether correct cluster was selected.')
    end
elseif strcmp(task_name, 'oddeven')      && strcmp(ref, 'target') && strcmp(cl_name, 'Beta')
    subclst = 50;
    % sanity check
    if cluster.subcluster(subclst).clsize(cluster.subcluster(subclst).largest) ~= 31 || ...
            cluster.subcluster(subclst).tmass(cluster.subcluster(subclst).largest) > -171.51 || ...
            cluster.subcluster(subclst).tmass(cluster.subcluster(subclst).largest) < -171.52
        warning('Check whether correct cluster was selected.')
    end
elseif strcmp(task_name, 'numberletter') && strcmp(ref, 'target') && strcmp(cl_name, 'Beta')
    subclst = 22;
    % sanity check
    if cluster.subcluster(subclst).clsize(cluster.subcluster(subclst).largest) ~= 23 || ...
            cluster.subcluster(subclst).tmass(cluster.subcluster(subclst).largest) > -89.82 || ...
            cluster.subcluster(subclst).tmass(cluster.subcluster(subclst).largest) < -89.83
        warning('Check whether correct cluster was selected.')
    end
end

% extract mean power values for cluster
voxels = cluster.subcluster(subclst).clusters == cluster.subcluster(subclst).largest;

% calculate mean power within subcluster
mean_power = nan(size(all_eegpower,1), size(all_eegpower,5));
for n = 1:size(all_eegpower,1)
    for ci = 1:size(all_eegpower,5)
        tmp = squeeze(all_eegpower(n,:,:,:,ci));
        mean_power(n,ci) = mean(tmp(voxels));
    end
end


%% ---------- Export ------------------------------------------------------

% export mean power values and subject list
mean_power = array2table(mean_power);
mean_power.Properties.VariableNames = {...
    [ref '_' cl_name '_' task_name '_switch_all'], ...
    [ref '_' cl_name '_' task_name '_repeat_all'], ...
    [ref '_' cl_name '_' task_name '_switch_odd'], ...
    [ref '_' cl_name '_' task_name '_switch_even'], ...  
    [ref '_' cl_name '_' task_name '_repeat_odd'], ...
    [ref '_' cl_name '_' task_name '_repeat_even'] }; ...
    % IMPORTANT NOTE: In matlab script "MUA_TF_03_wavelets.m" and its
    % output files, condition labels of 'switch_even' and 'repeat_odd' 
    % have been mixed up. In this script and its output files condition 
    % labels are correct.
mean_power = addvars(mean_power, char({allData.name}'), ...
    'Before', [ref '_' cl_name '_' task_name '_switch_all']);
writetable(mean_power, fullfile(PATH_IN, FILE_OUT_POWR))

