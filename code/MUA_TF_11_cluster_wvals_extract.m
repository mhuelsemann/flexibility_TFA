%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Title:         mass-univariate analysis with cluster-based permutation  %
%                tests of time-frequency data                             %
% Subtitle:      define core region of a cluster and extract stats        %
% Author:        Mareike J. Hülsemann, University of Mainz                %
% Date:          11/2023                                                  %
% Last Modified: 09/2025                                                  %
% Note:          Exp23, Switching Tasks                                   %
%                cluster pixels connected in space, frequency, and time   %
%                --> clustering for 3D-data with electrodes in first      %
%                dimension                                                %
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

measure   = 'wval';  % which measure is used for defining sub-clusters
OPTION    = '1SD'; logarithmize = true; % for core region
p         = '.005';  % significance level for parametric stats
alpha     = .05;     % for defining significance of a given cluster

% define PATHs
[filepath,~,~] = fileparts(matlab.desktop.editor.getActiveFilename);
parts          = regexp(filepath,filesep,'split');
PATH_MAIN      = fullfile(parts{1:end-1});

PATH_IN   = fullfile(PATH_MAIN, 'derivatives', task_name, ref, 'tf_stats');
PATH_OUT  = PATH_IN;

% define FILEs
FILE_IN_CLUST_ORIG = ['Exp23_' task_name '_tf_power_cluster_p' p(2:end) '.mat'];
FILE_OUT_PRINT     = ['Exp23_' task_name '_tf_power_cluster_p' p(2:end) '_wvals_' OPTION '.txt'];


%% ---------- Definition of Time-Frequency Cluster ------------------------

% load original cluster data
load(fullfile(PATH_IN, FILE_IN_CLUST_ORIG))

% define neighbouring elcs
max_elec_dist = .7;


% open text file to save results for publication
fileID = fopen(fullfile(PATH_IN, FILE_OUT_PRINT), 'w');

% Clustering for wvals
for i=1:cluster.amount_cluster
    
    % test whether current cluster is significant
    if cluster.pvals.tmass_overall(i) < alpha && ...
            cluster.pvals.clsize_overall(i) < alpha

        %% define min_wval

        % extract wvals
        [contr_elecs(:,1), contr_elecs(:,2), contr_elecs(:,3)] = ...
            ind2sub(size(cluster.clusters), find(cluster.clusters == i));
        wvals = nan(size(contr_elecs,1), 1);
        for j=1:size(contr_elecs,1)
            wvals(j,1) = cluster.TTEST(contr_elecs(j,1), contr_elecs(j,2), contr_elecs(j,3), 12);
        end
        
        % define cut-off wval for sub-cluster        
        if strcmp(OPTION((end-1):end), 'SD')  % OPTION 1: larger than MEAN + x SD
            sd2use = str2double(OPTION(1:end-2));
            if logarithmize
                min_wval = exp(mean(log(wvals)) + sd2use*std(log(wvals)));
            else
                min_wval = mean(wvals) + sd2use*std(wvals);
            end
        elseif strcmp(OPTION, 'upperthird')  % OPTION 2: > X % of wvals
            cutoff_idx = round(numel(wvals)*percentage);
            sorted_wvals = sort(wvals);
            min_wval = sorted_wvals(cutoff_idx+1);
        else  % OPTION 0: defining according to effect size conventions
        end

        clear contr_elecs wvals


        %% define sub-clusters

        % define sub-cluster according to w-val
        voxels = (cluster.clusters == i & cluster.TTEST(:,:,:,12) > min_wval);

        % test whether cluster splits up = defining new clusters

        % initialise cluster results
        tf_clust = zeros(size(voxels));

        % define clusters within electrodes
        cli=0;
        for ei=1:numel(cluster.info.chanlocs)
            for fi=1:numel(cluster.info.frex)
                for ti=1:numel(cluster.info.trial_time)

                    if fi==1 && ti==1 && voxels(ei,fi,ti) == 1
                        % begin of a cluster
                        cli=cli+1;
                        tf_clust(ei,fi,ti) = cli;
                    elseif fi~=1 && ti==1
                        if voxels(ei,fi,ti) == 1 && voxels(ei,fi-1,ti) == 0
                            % begin of a cluster
                            cli=cli+1;
                            tf_clust(ei,fi,ti) = cli;
                        elseif voxels(ei,fi,ti) == 1 && voxels(ei,fi-1,ti) == 1
                            % continuation of cluster
                            tf_clust(ei,fi,ti) = tf_clust(ei,fi-1,ti);
                        end
                    elseif fi==1 && ti~=1
                        if voxels(ei,fi,ti) == 1 && voxels(ei,fi,ti-1) == 0
                            % begin of a cluster
                            cli=cli+1;
                            tf_clust(ei,fi,ti) = cli;
                        elseif voxels(ei,fi,ti) == 1 && voxels(ei,fi,ti-1) == 1
                            % continuation of cluster
                            tf_clust(ei,fi,ti) = tf_clust(ei,fi,ti-1);
                        end
                    else % fi!=1 && ti!=1
                        if voxels(ei,fi,ti) == 1 && voxels(ei,fi,ti-1) == 0 && voxels(ei,fi-1,ti) == 0
                            % begin of a cluster
                            cli=cli+1;
                            tf_clust(ei,fi,ti) = cli;
                        elseif voxels(ei,fi,ti) == 1 && voxels(ei,fi,ti-1) == 1
                            % continuation of cluster
                            if voxels(ei,fi-1,ti) == 1
                                tf_clust(ei,fi,ti) = tf_clust(ei,fi-1,ti);
                                % rename previous cluster
                                tf_clust(tf_clust == tf_clust(ei,fi,ti-1)) = tf_clust(ei,fi-1,ti);
                            else
                                % continuation of cluster
                                tf_clust(ei,fi,ti) = tf_clust(ei,fi,ti-1);
                            end
                        elseif voxels(ei,fi,ti) == 1 && voxels(ei,fi-1,ti) == 1
                            % continuation of cluster
                            tf_clust(ei,fi,ti) = tf_clust(ei,fi-1,ti);
                        end
                    end
                end
            end
        end


        % define clusters across electrodes
        for ei=1:numel(cluster.info.chanlocs) % loop over all electrodes

            % get seed clusters
            tf_clust_seed = squeeze(tf_clust(ei,:,:));

            new_clusters = unique(tf_clust_seed);
            new_clusters(new_clusters == 0) = [];

            % loop over all elec_clusters of the current electrode
            for ci=1:numel(new_clusters)
                % define pixel idxs of current seed cluster
                idx_seed = find(tf_clust_seed == new_clusters(ci));

                % Are there neighbours with overlapping clusters in time?
                for ni=1:numel(cluster.info.chanlocs)
                    % find potential neighbours
                    if ei==ni % do nothing, as this is the same electrode
                        % verify, whether electrode is direct neighbour
                    elseif sqrt((cluster.info.chanlocs(ei).X - cluster.info.chanlocs(ni).X)^2 + ...
                            (cluster.info.chanlocs(ei).Y - cluster.info.chanlocs(ni).Y)^2 + ...
                            (cluster.info.chanlocs(ei).Z - cluster.info.chanlocs(ni).Z)^2) ...
                            < max_elec_dist

                        % get neighbours clusters
                        tf_clust_nb = squeeze(tf_clust(ni,:,:));
                        % define pixel idxs of all neighbour clusters
                        idx_nb = find(tf_clust_nb ~= 0);
                        % Any overlapping clusters?
                        overlap = intersect(idx_seed, idx_nb);

                        if ~isempty(overlap)
                            % names of overlapping clusters
                            cl_names = unique(tf_clust_nb(overlap));
                            % define minimal cluster number
                            min_clust_nr = min([cl_names; new_clusters(ci)]);

                            % rename overlapping clusters to minimal cluster number
                            tf_clust_seed(tf_clust_seed == new_clusters(ci)) = min_clust_nr;
                            tf_clust_nb(ismember(tf_clust_nb, cl_names)) = min_clust_nr;

                            % re-write cluster numbers in original data
                            tf_clust(ei,:,:) = tf_clust_seed;
                            tf_clust(ni,:,:) = tf_clust_nb;
                        end
                    end
                end
            end
        end

        % change cluster numbers to 1-by-1 increasing cluster numbers
        clusternames = unique(tf_clust);
        clusternames(clusternames==0) = [];
        for ci=1:numel(clusternames)
            tf_clust(tf_clust==clusternames(ci)) = ci;
        end

        % define all clusters
        all_clusters = unique(tf_clust);
        all_clusters(all_clusters==0) = [];
        
        % define cluster sizes
        tvals  = squeeze(cluster.TTEST(:,:,:,1));
        tmass  = nan(size(all_clusters));
        clsize = nan(size(all_clusters));
        for ci=1:numel(all_clusters) % loop through all clusters
            tmass(ci)  =   sum(tvals(tf_clust==all_clusters(ci)));
            clsize(ci) = numel(tvals(tf_clust==all_clusters(ci)));
        end

        % define largest sub-cluster and save results
        if find(abs(tmass)==max(abs(tmass))) == find(clsize==max(clsize))
            % unambigious cluster found
            cluster.subcluster(i).clusters = tf_clust;
            cluster.subcluster(i).tmass    = tmass;
            cluster.subcluster(i).clsize   = clsize;
            cluster.subcluster(i).largest  = find(clsize==max(clsize));
            cluster.subcluster(i).info.measure = measure;
            cluster.subcluster(i).info.option  = OPTION;
            cluster.subcluster(i).info.log     = logarithmize;
            cluster.subcluster(i).info.cutoff  = min_wval;
        else 
            % largest sub-cluster not identical to largest tmass-sub-cluster
            warning("Largest sub-cluster not identical to largest tmass-sub-cluster")
        end


        %% print out largst sub-cluster results

        subcluster = tf_clust == find(clsize==max(clsize));

        % extract stats for relevant cluster
        [contr_elecs(:,1), contr_elecs(:,2), contr_elecs(:,3)] = ...
            ind2sub(size(subcluster), find(subcluster));
        tvals = nan(size(contr_elecs,1), 1); 
        wvals = nan(size(contr_elecs,1), 1);
        for j=1:size(contr_elecs,1)
            tvals(j,1) = cluster.TTEST(contr_elecs(j,1), contr_elecs(j,2), contr_elecs(j,3),  1);
            wvals(j,1) = cluster.TTEST(contr_elecs(j,1), contr_elecs(j,2), contr_elecs(j,3), 12); % Omega^2_partial
        end
        unique_elecs = unique(contr_elecs(:,1));
        time_lb = cluster.info.trial_time(min(contr_elecs(:,3)));
        time_ub = cluster.info.trial_time(max(contr_elecs(:,3)));        
        frex_lb = cluster.info.frex(min(contr_elecs(:,2)));
        frex_ub = cluster.info.frex(max(contr_elecs(:,2)));
        max_elec = find(abs(tvals)==max(abs(tvals)));
                CI_lo = mean(wvals) - 1.96 * std(wvals)/sqrt(numel(wvals));
        CI_up = mean(wvals) + 1.96 * std(wvals)/sqrt(numel(wvals));
        
        % write cluster stats to file
        fprintf(fileID, 'Cluster %d is significant.\r\n',i);

        fprintf(fileID, '%s  &  %s,p=%s,w=%0.3f,%s (%0.2f,%d) &  %0.1f--%0.1f  &  %d--%d  & %0.2f  &  %d  &  %0.3f $\\pm$ %0.3f  & ', ...
            task_name, ref, p, min_wval, OPTION, ...
            cluster.tmass(i), cluster.clsize(i), ...
            frex_lb, frex_ub, time_lb, time_ub, ...
            cluster.subcluster(i).tmass(cluster.subcluster(i).largest), ...
            cluster.subcluster(i).clsize(cluster.subcluster(i).largest), ...
            mean(wvals), std(wvals));
        for j=1:numel(unique_elecs)
            fprintf(fileID, '%s ', cluster.info.chanlocs(unique_elecs(j)).labels);
        end
        fprintf(fileID, '\\\\\r\n\r\n');

        fprintf(fileID, 'The electrode with the largest effect is %s at %0.2f Hz and %d ms, t(%d) = %0.2f, p = %0.3f, $\\omega_p^2$ = %0.3f.\r\n', ...
            cluster.info.chanlocs(contr_elecs(max_elec,1)).labels, ...
            cluster.info.frex(contr_elecs(max_elec,2)), ... % frequency
            cluster.info.trial_time(contr_elecs(max_elec,3)), ... % timewindow
            cluster.TTEST(contr_elecs(max_elec,1), contr_elecs(max_elec,2), contr_elecs(max_elec,3),  2), ... % dfN
            cluster.TTEST(contr_elecs(max_elec,1), contr_elecs(max_elec,2), contr_elecs(max_elec,3),  1), ... % t
            cluster.TTEST(contr_elecs(max_elec,1), contr_elecs(max_elec,2), contr_elecs(max_elec,3),  3), ... % pval
            cluster.TTEST(contr_elecs(max_elec,1), contr_elecs(max_elec,2), contr_elecs(max_elec,3), 12));  ... % wval
        

        fprintf(fileID, '\r\n\r\n');

        clear contr_elecs tvals wvals max_elec

    end
end
fclose(fileID);


%% Save cluster results

% save
save(fullfile(PATH_OUT, FILE_IN_CLUST_ORIG),'cluster')


