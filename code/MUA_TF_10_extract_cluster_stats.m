%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Title:         mass-univariate analysis with cluster-based permutation  %
%                tests of time-frequency data                             %
% Subtitle:      extraction of cluster statistics                         %
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
task_name = 'numberletter';
ref       = 'target';

p         = '.005'; % significance level for parametric stats; which voxel is above a threshold to belong to a cluster
nperm     = 1000; % amount of permutations
alpha     = .05; % for defining significance of a given cluster

% define PATHs
[filepath,~,~] = fileparts(matlab.desktop.editor.getActiveFilename);
parts          = regexp(filepath,filesep,'split');
PATH_MAIN      = fullfile(parts{1:end-1});

PATH_IN   = fullfile(PATH_MAIN, 'derivatives', task_name, ref, 'tf_stats');
PATH_OUT  = PATH_IN;

% define FILEs
FILE_IN_TTEST      = ['Exp23_' task_name '_tf_power_TTEST.mat'];
FILE_IN_CLUST_ORIG = ['Exp23_' task_name '_tf_power_cluster_p' p(2:end) '.mat'];
FILE_IN_CLUST_PERM = ['Exp23_' task_name '_tf_power_cluster_p' p(2:end) '_perm.mat'];
FILE_OUT_PRINT     = ['Exp23_' task_name '_tf_power_cluster_p' p(2:end) '_stats_'];


%% ---------- Load data ---------------------------------------------------

% load original and permuted cluster data
load(fullfile(PATH_IN, FILE_IN_CLUST_ORIG))
load(fullfile(PATH_IN, FILE_IN_CLUST_PERM))

% load original parametric stats
load(fullfile(PATH_IN, FILE_IN_TTEST))


%% ---------- Calculate cluster stats -------------------------------------

% initialise
pvals_varnames = {'clsize_overall', 'tmass_overall', 'clsize_pos', 'tmass_pos', 'clsize_neg', 'tmass_neg'}; 
pvals     = nan(cluster.amount_cluster, numel(pvals_varnames));

% calculate p-val
for i=1:cluster.amount_cluster
    pvals(i,1) = sum( cluster_permuted.max_clsize_overall > cluster.clsize(i)      ) / nperm;
    pvals(i,2) = sum( cluster_permuted.max_tmass_overall  > abs(cluster.tmass(i))  ) / nperm;
    
    if cluster.tmass(i) >= 0
        pvals(i,3) = sum( cluster_permuted.max_clsize_pos > cluster.clsize(i) ) / nperm;
        pvals(i,4) = sum( cluster_permuted.max_tmass_pos  > cluster.tmass(i)  ) / nperm;
    else
        pvals(i,5) = sum( cluster_permuted.max_clsize_neg > cluster.clsize(i) ) / nperm;
        pvals(i,6) = sum( cluster_permuted.max_tmass_neg  < cluster.tmass(i)  ) / nperm;
    end
    
end

% organise pvals in concise table
pvals = array2table(pvals);
pvals.Properties.VariableNames = pvals_varnames;
cluster_sign(cluster.tmass >= 0) = 1;
cluster_sign(cluster.tmass <  0) = -1;
pvals = addvars(pvals, [1:cluster.amount_cluster]', cluster_sign', ...
    'Before', 'clsize_overall', ...
    'NewVariableNames', {'id', 'tmass_sign'});

% add pvals to cluster structure
cluster.pvals = pvals;

% save
save(fullfile(PATH_IN, FILE_IN_CLUST_ORIG),'cluster')


%% ---------- print out stats for each cluster (overall) ------------------

% open text file to save results for publication
fileID = fopen(fullfile(PATH_IN, [FILE_OUT_PRINT 'overall.txt']), 'w');
for i=1:cluster.amount_cluster
    
    % test whether current cluster is significant
    if cluster.pvals.tmass_overall(i) >= alpha && ...
            cluster.pvals.clsize_overall(i) >= alpha
        fprintf(fileID, 'Cluster %d is not significant.\r\n\r\n',i);
    else
        fprintf(fileID, 'Cluster %d is significant.\r\n',i);
        fprintf(fileID, 't-mass: %0.2f, p=%0.3f;\r\ncluster size: %d, p=%0.3f.\r\n', ...
            cluster.tmass(i), cluster.pvals.tmass_overall(i), ...
            cluster.clsize(i), cluster.pvals.clsize_overall(i));
        fprintf(fileID, 'The following electrodes contribute to the cluster: ');
        
        [contr_elecs(:,1), contr_elecs(:,2), contr_elecs(:,3)] = ...
            ind2sub(size(cluster.clusters), find(cluster.clusters == i));        
        tvals = nan(size(contr_elecs,1), 1); 
        wvals = nan(size(contr_elecs,1), 1);
        for j=1:size(contr_elecs,1)
            tvals(j,1) = cluster.TTEST(contr_elecs(j,1), contr_elecs(j,2), contr_elecs(j,3),  1);
            wvals(j,1) = cluster.TTEST(contr_elecs(j,1), contr_elecs(j,2), contr_elecs(j,3), 12); % Omega^2_partial
%             fprintf(fileID, '%s, %0.2f Hz, %d ms: t(%d) = %0.2f, p = %0.3f, omega_p^2 = %0.3f (r = %0.2f)\r\n', ...
%                 cluster.info.chanlocs(contr_elecs(j,1)).labels, ... % elec
%                 cluster.info.frex(contr_elecs(j,2)), ... % freq
%                 cluster.info.trial_time(contr_elecs(j,3)), ... % timewindow
%                 cluster.TTEST(contr_elecs(j,1), contr_elecs(j,2), contr_elecs(j,3),  2), ... % dfN
%                 cluster.TTEST(contr_elecs(j,1), contr_elecs(j,2), contr_elecs(j,3),  1), ... % t
%                 cluster.TTEST(contr_elecs(j,1), contr_elecs(j,2), contr_elecs(j,3),  3), ... % pval
%                 cluster.TTEST(contr_elecs(j,1), contr_elecs(j,2), contr_elecs(j,3), 12), ... % wval
%                 cluster.TTEST(contr_elecs(j,1), contr_elecs(j,2), contr_elecs(j,3),  4));    % r
        end
        unique_elecs = unique(contr_elecs(:,1));
        for j=1:numel(unique_elecs)
            fprintf(fileID, '%s ', cluster.info.chanlocs(unique_elecs(j)).labels);
        end
        fprintf(fileID, '\r\n');
        
        time_lb = trial_time(min(contr_elecs(:,3)));
        time_ub = trial_time(max(contr_elecs(:,3)));        
        frex_lb = frex(min(contr_elecs(:,2)));
        frex_ub = frex(max(contr_elecs(:,2)));
        fprintf(fileID, 'We detected a cluster of significant differences between conditions %d-%d ms after %s onset for frequencies between %0.1f-%0.1f Hz (t-mass = %0.2f, p = %0.3f; cluster size = %d, p = %0.3f).\r\n', ...
            time_lb, time_ub, ref, frex_lb, frex_ub, ...
            cluster.tmass(i), cluster.pvals.tmass_overall(i), ...
            cluster.clsize(i), cluster.pvals.clsize_overall(i));
                
        max_elec = find(abs(tvals)==max(abs(tvals)));
        fprintf(fileID, 'The electrode with the largest effect is %s at %0.2f Hz and %d ms, t(%d) = %0.2f, p = %0.3f, omega_p^2 = %0.3f.\r\n', ...
            cluster.info.chanlocs(contr_elecs(max_elec,1)).labels, ...
            cluster.info.frex(contr_elecs(max_elec,2)), ... % frequency
            cluster.info.trial_time(contr_elecs(max_elec,3)), ... % timewindow
            cluster.TTEST(contr_elecs(max_elec,1), contr_elecs(max_elec,2), contr_elecs(max_elec,3),  2), ... % dfN
            cluster.TTEST(contr_elecs(max_elec,1), contr_elecs(max_elec,2), contr_elecs(max_elec,3),  1), ... % t
            cluster.TTEST(contr_elecs(max_elec,1), contr_elecs(max_elec,2), contr_elecs(max_elec,3),  3), ... % pval
            cluster.TTEST(contr_elecs(max_elec,1), contr_elecs(max_elec,2), contr_elecs(max_elec,3), 12));  ... % wval
        
        CI_lo = mean(wvals) - 1.96 * std(wvals)/sqrt(numel(wvals));
        CI_up = mean(wvals) + 1.96 * std(wvals)/sqrt(numel(wvals));
        fprintf(fileID, ...
            '(t-mass  = %0.2f, p = %.3f; cluster size = %d, p = %.3f; mean omega_p^2 = %.3f, 95%% CI [%.3f, %.3f])\r\n', ...
            cluster.tmass(i), cluster.pvals.tmass_overall(i), ...
            cluster.clsize(i), cluster.pvals.clsize_overall(i), ...
            mean(wvals), CI_lo, CI_up);
        fprintf(fileID, ...
            'Effect sizes for the contributing electrodes ranged between %0.3f < omega_p^2 < %0.3f, with an average (plusminus SD) of omega_p^2 = %0.3f plusminus %0.3f.\r\n\r\n', ...
            min(wvals), max(wvals), mean(wvals), std(wvals));
        
        fprintf(fileID, '%s  &  p=%s  &  %0.1f--%0.1f  &  %d--%d  & %0.2f  &  %0.3f  &  %d  &  %0.3f  &  %0.3f $\\pm$ %0.3f  & ', ...
            task_name, p, frex_lb, frex_ub, time_lb, time_ub, ...
            cluster.tmass(i), cluster.pvals.tmass_overall(i), ...
            cluster.clsize(i), cluster.pvals.clsize_overall(i), ...
            mean(wvals), std(wvals));
        for j=1:numel(unique_elecs)
            fprintf(fileID, '%s ', cluster.info.chanlocs(unique_elecs(j)).labels);
        end
        fprintf(fileID, '\\\\\r\n\r\n');

        clear contr_elecs tvals wvals max_elec
    end
end
fclose(fileID);