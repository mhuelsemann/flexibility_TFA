%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Title:         mass-univariate analysis with cluster-based permutation  %
%                tests of time-frequency data                             %
% Subtitle:      cluster empirical data                                   %
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

note      = 'same amount of trials';
p         = '.005';


% define PATHs
[filepath,~,~] = fileparts(matlab.desktop.editor.getActiveFilename);
parts          = regexp(filepath,filesep,'split');
PATH_MAIN      = fullfile(parts{1:end-1});

PATH_IN   = fullfile(PATH_MAIN, 'derivatives', task_name, ref, 'tf_stats');
PATH_OUT  = PATH_IN;

% define FILEs
FILE_IN  = ['Exp23_' task_name '_tf_power_TTEST.mat'];
FILE_OUT = ['Exp23_' task_name '_tf_power_cluster_p' p(2:end) '.mat'];



%% ---------- Definition of Time-Frequency Cluster ------------------------

% load empirical t-vals
load(fullfile(PATH_IN, FILE_IN))
tvals =        squeeze(TTEST(:,:,:,1));
df    = unique(squeeze(TTEST(:,:,:,2)));

% define critical t-val
tails = 'twosided';
alpha = str2double(p);
if strcmp(tails,'twosided'), alpha = alpha/2; end
tcrit = tinv(1-alpha,df); % critical t-value

% define neighbouring elcs
max_elec_dist = .7;


%% Clustering for positive t-vals

% define (positive) significant t-vals
tvals_psig = zeros(size(tvals));
tvals_psig(tvals > tcrit) = 1;

% initialise (positive) cluster results
tf_pclust = zeros(size(tvals_psig));

% define (positive) clusters within electrodes
cli=0;
for ei=1:numel(chanlocs)
    for fi=1:size(tvals_psig,2)
        for ti=1:size(tvals_psig,3)
            
            if fi==1 && ti==1 && tvals_psig(ei,fi,ti) == 1
                % begin of a cluster
                cli=cli+1;
                tf_pclust(ei,fi,ti) = cli;
            elseif fi~=1 && ti==1
                if tvals_psig(ei,fi,ti) == 1 && tvals_psig(ei,fi-1,ti) == 0
                    % begin of a cluster
                    cli=cli+1;
                    tf_pclust(ei,fi,ti) = cli;
                elseif tvals_psig(ei,fi,ti) == 1 && tvals_psig(ei,fi-1,ti) == 1
                    % continuation of cluster
                    tf_pclust(ei,fi,ti) = tf_pclust(ei,fi-1,ti);
                end
            elseif fi==1 && ti~=1
                if tvals_psig(ei,fi,ti) == 1 && tvals_psig(ei,fi,ti-1) == 0
                    % begin of a cluster
                    cli=cli+1;
                    tf_pclust(ei,fi,ti) = cli;
                elseif tvals_psig(ei,fi,ti) == 1 && tvals_psig(ei,fi,ti-1) == 1
                    % continuation of cluster
                    tf_pclust(ei,fi,ti) = tf_pclust(ei,fi,ti-1);
                end
            else % fi!=1 && ti!=1
                if tvals_psig(ei,fi,ti) == 1 && tvals_psig(ei,fi,ti-1) == 0 && tvals_psig(ei,fi-1,ti) == 0
                    % begin of a cluster
                    cli=cli+1;
                    tf_pclust(ei,fi,ti) = cli;
                elseif tvals_psig(ei,fi,ti) == 1 && tvals_psig(ei,fi,ti-1) == 1
                    % continuation of cluster
                    if tvals_psig(ei,fi-1,ti) == 1
                        tf_pclust(ei,fi,ti) = tf_pclust(ei,fi-1,ti);
                        % rename previous cluster
                        tf_pclust(tf_pclust == tf_pclust(ei,fi,ti-1)) = tf_pclust(ei,fi-1,ti);
                    else
                        % continuation of cluster
                        tf_pclust(ei,fi,ti) = tf_pclust(ei,fi,ti-1);
                    end
                elseif tvals_psig(ei,fi,ti) == 1 && tvals_psig(ei,fi-1,ti) == 1
                    % continuation of cluster
                    tf_pclust(ei,fi,ti) = tf_pclust(ei,fi-1,ti);
                end
            end
        end
    end
end

% define (positive) clusters across electrodes
for ei=1:numel(chanlocs) % loop over all electrodes
    
    % get seed clusters
    tf_clust_seed = squeeze(tf_pclust(ei,:,:));
    
    clusters = unique(tf_clust_seed);
    clusters(clusters == 0) = [];
    
    % loop over all elec_clusters of the current electrode
    for ci=1:numel(clusters)
        % define pixel idxs of current seed cluster
        idx_seed = find(tf_clust_seed == clusters(ci));
        
        % Are there neighbours with overlapping clusters in time?
        for ni=1:numel(chanlocs)
            % find potential neighbours
            if ei==ni % do nothing, as this is the same electrode
            % verify, whether electrode is direct neighbour
            elseif sqrt((chanlocs(ei).X - chanlocs(ni).X)^2 + ...
                    (chanlocs(ei).Y - chanlocs(ni).Y)^2 + ...
                    (chanlocs(ei).Z - chanlocs(ni).Z)^2) < max_elec_dist
                
                % get neighbours clusters
                tf_clust_nb = squeeze(tf_pclust(ni,:,:));
                % define pixel idxs of all neighbour clusters
                idx_nb = find(tf_clust_nb ~= 0);
                % Any overlapping clusters?
                overlap = intersect(idx_seed, idx_nb);
                
                if ~isempty(overlap)
                    % names of overlapping clusters
                    cl_names = unique(tf_clust_nb(overlap));
                    % define minimal cluster number
                    min_clust_nr = min([cl_names; clusters(ci)]);
                    
                    % rename overlapping clusters to minimal cluster number
                    tf_clust_seed(tf_clust_seed == clusters(ci)) = min_clust_nr;
                    tf_clust_nb(ismember(tf_clust_nb, cl_names)) = min_clust_nr;
                    
                    % re-write cluster numbers in original data
                    tf_pclust(ei,:,:) = tf_clust_seed;
                    tf_pclust(ni,:,:) = tf_clust_nb;
                end
            end
        end
    end
end


%% Clustering for for negative t-vals

% define (negative) significant t-vals
tvals_nsig = zeros(size(tvals));
tvals_nsig(tvals < (-tcrit)) = 1;

% initialise (negative) cluster results
tf_nclust = zeros(size(tvals_nsig));

% define (negative) clusters within electrodes
cli=max(tf_pclust(:))+1;
for ei=1:numel(chanlocs)
    for fi=1:size(tvals_nsig,2)
        for ti=1:size(tvals_nsig,3)
            
            if fi==1 && ti==1 && tvals_nsig(ei,fi,ti) == 1
                % begin of a cluster
                cli=cli+1;
                tf_nclust(ei,fi,ti) = cli;
            elseif fi~=1 && ti==1
                if tvals_nsig(ei,fi,ti) == 1 && tvals_nsig(ei,fi-1,ti) == 0
                    % begin of a cluster
                    cli=cli+1;
                    tf_nclust(ei,fi,ti) = cli;
                elseif tvals_nsig(ei,fi,ti) == 1 && tvals_nsig(ei,fi-1,ti) == 1
                    % continuation of cluster
                    tf_nclust(ei,fi,ti) = tf_nclust(ei,fi-1,ti);
                end
            elseif fi==1 && ti~=1
                if tvals_nsig(ei,fi,ti) == 1 && tvals_nsig(ei,fi,ti-1) == 0
                    % begin of a cluster
                    cli=cli+1;
                    tf_nclust(ei,fi,ti) = cli;
                elseif tvals_nsig(ei,fi,ti) == 1 && tvals_nsig(ei,fi,ti-1) == 1
                    % continuation of cluster
                    tf_nclust(ei,fi,ti) = tf_nclust(ei,fi,ti-1);
                end
            else % fi!=1 && ti!=1
                if tvals_nsig(ei,fi,ti) == 1 && tvals_nsig(ei,fi,ti-1) == 0 && tvals_nsig(ei,fi-1,ti) == 0
                    % begin of a cluster
                    cli=cli+1;
                    tf_nclust(ei,fi,ti) = cli;
                elseif tvals_nsig(ei,fi,ti) == 1 && tvals_nsig(ei,fi,ti-1) == 1
                    % continuation of cluster
                    if tvals_nsig(ei,fi-1,ti) == 1
                        tf_nclust(ei,fi,ti) = tf_nclust(ei,fi-1,ti);
                        % rename previous cluster
                        tf_nclust(tf_nclust == tf_nclust(ei,fi,ti-1)) = tf_nclust(ei,fi-1,ti);
                    else
                        % continuation of cluster
                        tf_nclust(ei,fi,ti) = tf_nclust(ei,fi,ti-1);
                    end
                elseif tvals_nsig(ei,fi,ti) == 1 && tvals_nsig(ei,fi-1,ti) == 1
                    % continuation of cluster
                    tf_nclust(ei,fi,ti) = tf_nclust(ei,fi-1,ti);
                end
            end
        end
    end
end

% define (negative) clusters across electrodes
for ei=1:numel(chanlocs) % loop over all electrodes
    
    % get seed clusters
    tf_clust_seed = squeeze(tf_nclust(ei,:,:));
    
    clusters = unique(tf_clust_seed);
    clusters(clusters == 0) = [];
    
    % loop over all elec_clusters of the current electrode
    for ci=1:numel(clusters)
        % define pixel idxs of current seed cluster
        idx_seed = find(tf_clust_seed == clusters(ci));
        
        % Are there neighbours with overlapping clusters in time?
        for ni=1:numel(chanlocs)
            % find potential neighbours
            if ei==ni % do nothing, as this is the same electrode
            % verify, whether electrode is direct neighbour
            elseif sqrt((chanlocs(ei).X - chanlocs(ni).X)^2 + ...
                    (chanlocs(ei).Y - chanlocs(ni).Y)^2 + ...
                    (chanlocs(ei).Z - chanlocs(ni).Z)^2) < max_elec_dist
                
                % get neighbours clusters
                tf_clust_nb = squeeze(tf_nclust(ni,:,:));
                % define pixel idxs of all neighbour clusters
                idx_nb = find(tf_clust_nb ~= 0);
                % Any overlapping clusters?
                overlap = intersect(idx_seed, idx_nb);
                
                if ~isempty(overlap)
                    % names of overlapping clusters
                    cl_names = unique(tf_clust_nb(overlap));
                    % define minimal cluster number
                    min_clust_nr = min([cl_names; clusters(ci)]);
                    
                    % rename overlapping clusters to minimal cluster number
                    tf_clust_seed(tf_clust_seed == clusters(ci)) = min_clust_nr;
                    tf_clust_nb(ismember(tf_clust_nb, cl_names)) = min_clust_nr;
                    
                    % re-write cluster numbers in original data
                    tf_nclust(ei,:,:) = tf_clust_seed;
                    tf_nclust(ni,:,:) = tf_clust_nb;
                end
            end
        end
    end
end


%% Summarize cluster results

% merge positive and negative cluster
tf_clust = tf_pclust + tf_nclust;

% change cluster numbers to 1-by-1 increasing cluster numbers
clusternames = unique(tf_clust);
clusternames(clusternames==0) = [];
for i=1:numel(clusternames)
    tf_clust(tf_clust==clusternames(i)) = i;
end

% define all clusters
all_clusters = unique(tf_clust);
all_clusters(all_clusters==0) = [];
    
% define cluster sizes
tmass  = nan(size(all_clusters));
clsize = nan(size(all_clusters));
for ci=1:numel(all_clusters) % loop through all clusters
    tmass(ci)  =   sum(tvals(tf_clust==all_clusters(ci)));
    clsize(ci) = numel(tvals(tf_clust==all_clusters(ci)));
end


%% ---------- Save results ------------------------------------------------

% save all cluster information
cluster.amount_cluster = numel(all_clusters);
cluster.clusters       = tf_clust;
cluster.tmass          = tmass;
cluster.clsize         = clsize;

% save information about data structure
info.chanlocs     = chanlocs;
info.frex         = frex;
info.trial_time   = trial_time;
info.clust_struct = 'electrodes x frequencies x time';
info.TTEST_file   = fullfile(PATH_IN, FILE_IN);
info.alpha        = 2*alpha;
info.tails        = tails;
info.tcrit        = tcrit;
info.note         = note;
% add to cluster_result structure
cluster.info = info;

% add parametric stats to cluster
cluster.TTEST           = TTEST;
cluster.info.TTEST_vars = TTEST_vars;
cluster.info.allData    = allData;

% save
save(fullfile(PATH_OUT, FILE_OUT),'cluster')
