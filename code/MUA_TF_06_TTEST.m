%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Title:         mass-univariate analysis with cluster-based permutation  %
%                tests of time-frequency data                             %
% Subtitle:      t-tests for empirical data                               %
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
load(fullfile(PATH_IN, allData(1).name), 'chanlocs', 'frex', 'trial_time')
all_eegpower = nan(numel(allData),numel(chanlocs),numel(frex),numel(trial_time),2);
for n=1:numel(allData)
    load(fullfile(PATH_IN, allData(n).name), 'eegpower');
    all_eegpower(n,:,:,:,:) = eegpower(:,:,:,1:2);
end
clear eegpower


%% ---------- Execute t-Tests ---------------------------------------------

% initialise
TTEST_vars = {'T', 'df', 'p', 'r', ...
    'd_z','d_rm','g_rm','d_av','g_av','CL', ...
    'omega_within', 'omega_p_within'};
TTEST = nan(numel(chanlocs), numel(frex), numel(trial_time), numel(TTEST_vars));

for ci=1:numel(chanlocs)
    tmp = squeeze(all_eegpower(:,ci,:,:,:));
    tmp_TTEST = nan(numel(frex), numel(trial_time), numel(TTEST_vars));
    
    for fi=1:numel(frex)
        for ti=1:numel(trial_time)
            
            % get data
            data1 = squeeze(tmp(:,fi,ti,1));
            data2 = squeeze(tmp(:,fi,ti,2));
            data  = [data1, data2];
            
            % run t-tests
            [H,P,CI,STATS] = ttest(data(:,1), data(:,2));
            Tval = STATS.tstat;
            df   = STATS.df;
            N    = size(data,1);
            nobs = 2*N;
            
            % calc difference scores
            diff    = data(:,1) - data(:,2);
            diff_mn = abs(mean(diff));
            diff_sd = std(diff); % = sqrt( std(data(:,1))^2 + std(data(:,2))^2 - (2*r*std(data(:,1))*std(data(:,2))) );
            
            % calculate d-effect sizes
            d_z     = diff_mn/diff_sd; % = abs(Tval)/sqrt(N);
            r       = corr(data); r = r(1,2);
            d_rm    = d_z * sqrt(2*(1-r)); % = d_rm_rsign;
            g_rm    = d_rm * (1-(3/(8*df-1))); % = g_rm_rsign;
            d_av    = diff_mn / ( (std(data(:,1))+std(data(:,2)))/2 );
            g_av    = d_av*(1-(3/(8*df-1)));
            CL      = cdf(makedist('Normal'), diff_mn/diff_sd); % from Lakens spreadsheet
                        
            % means, SS, df
            subjects_mean = mean(data, 2);
            condt_mean = mean(data, 1);
            total_mean = mean(condt_mean);
            
            SS_total   = sum(sum((data - total_mean).^2));
            SS_effect = sum((condt_mean - total_mean).^2 * N);
            SS_error = sum(sum((data - repmat(condt_mean,N,1) - [subjects_mean, subjects_mean] + total_mean).^2));
            SS_subjects   = sum((subjects_mean - total_mean).^2 * 2); % p = 2 for t-tests
            
            df_total   = 2*N - 1; % pn - 1
            df_effect = 1;       % p  - 1
            df_error = N - 1;   % (p-1)(n-1)
            df_subjects   = N - 1;   %  n - 1
            
            MS_total   = SS_total   / df_total;   % total
            MS_effect = SS_effect / df_effect; % effect
            MS_error = SS_error / df_error; % effect*subjects
            MS_subjects   = SS_subjects   / df_subjects;   % subjects
            
            F = MS_effect / MS_error;
            
            % calculate effect size omega^2
            omega_within = (df_effect * (MS_effect-MS_error)) / (SS_total + MS_subjects);
            fs   = ((Tval^2)-1) / nobs;    % only for within-subjects design
            omega_p_within = (fs)/(1+fs);  % " " "
                        
            col =     1; tmp_TTEST(fi,ti,col) = Tval;
            col = col+1; tmp_TTEST(fi,ti,col) = df;
            col = col+1; tmp_TTEST(fi,ti,col) = P;
            col = col+1; tmp_TTEST(fi,ti,col) = r;
            col = col+1; tmp_TTEST(fi,ti,col) = d_z;
            col = col+1; tmp_TTEST(fi,ti,col) = d_rm;
            col = col+1; tmp_TTEST(fi,ti,col) = g_rm;
            col = col+1; tmp_TTEST(fi,ti,col) = d_av;
            col = col+1; tmp_TTEST(fi,ti,col) = g_av;
            col = col+1; tmp_TTEST(fi,ti,col) = CL;
            col = col+1; tmp_TTEST(fi,ti,col) = omega_within;
            col = col+1; tmp_TTEST(fi,ti,col) = omega_p_within;
           
        end
    end

    TTEST(ci,:,:,:) = tmp_TTEST;
end


% save output
save(fullfile(PATH_OUT, [allData(1).name(1:6) task_name '_tf_power_TTEST.mat']), ...
    'all_eegpower', 'allData', 'TTEST', 'note', 'task_name', 'ref', ...
    'chanlocs', 'frex', 'trial_time', 'TTEST_vars');

