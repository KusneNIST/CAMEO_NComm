function [al_results] = run_active_learning_experiment_190618a2(experiment)

% This software was developed by employees of the National Institute of
% Standards and Technology (NIST), an agency of the Federal Government and
% is being made available as a public service. Pursuant to title 17 United
% States Code Section 105, works of NIST employees are not subject to
% copyright protection in the United States.  This software may be subject
% to foreign copyright.  Permission in the United States and in foreign
% countries, to the extent that NIST may hold copyright, to use, copy,
% modify, create derivative works, and distribute this software and its
% documentation without fee is hereby granted on a non-exclusive basis,
% provided that this notice and disclaimer of warranty appears in all
% copies.

% THE SOFTWARE IS PROVIDED 'AS IS' WITHOUT ANY WARRANTY OF ANY KIND, EITHER
% EXPRESSED, IMPLIED, OR STATUTORY, INCLUDING, BUT NOT LIMITED TO, ANY
% WARRANTY THAT THE SOFTWARE WILL CONFORM TO SPECIFICATIONS, ANY IMPLIED
% WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, AND
% FREEDOM FROM INFRINGEMENT, AND ANY WARRANTY THAT THE DOCUMENTATION WILL
% CONFORM TO THE SOFTWARE, OR ANY WARRANTY THAT THE SOFTWARE WILL BE ERROR
% FREE.  IN NO EVENT SHALL NIST BE LIABLE FOR ANY DAMAGES, INCLUDING, BUT
% NOT LIMITED TO, DIRECT, INDIRECT, SPECIAL OR CONSEQUENTIAL DAMAGES,
% ARISING OUT OF, RESULTING FROM, OR IN ANY WAY CONNECTED WITH THIS
% SOFTWARE, WHETHER OR NOT BASED UPON WARRANTY, CONTRACT, TORT, OR
% OTHERWISE, WHETHER OR NOT INJURY WAS SUSTAINED BY PERSONS OR PROPERTY OR
% OTHERWISE, AND WHETHER OR NOT LOSS WAS SUSTAINED FROM, OR AROSE OUT OF
% THE RESULTS OF, OR USE OF, THE SOFTWARE OR SERVICES PROVIDED HEREUNDER.

% A. Gilad Kusne, NIST, aaron.kusne@nist.gov, Release 8/01/2020
% If using this work for a publication, please cite:
% Kusne, A. G., et al. "On-the-fly Closed-loop Autonomous Materials
% Discovery via Bayesian Active Learning." arXiv preprint arXiv:2006.06141
% (2020).

% load the appropriate data if simulation or real experiment.
[X, T, C, XY, ~, experiment] = load_wafer_data(experiment);
experiment.point_position_to_index = [(1:size(XY,1))' XY]; %
experiment.performance.fmiw = [];
Cexp_XY = tern2cart(C);

% Select initial points to measure
coarse_indices = experiment.coarse_indices;

% initial measurements
Xa = X(coarse_indices, :);
Fa = experiment.functional_property(coarse_indices);

% composition data to graph
S = form_graph_171116a(C, XY, experiment);
S = use_prior_to_impact_S(experiment, S);

measured_indices_just_experimental = coarse_indices';
num_new = length(measured_indices_just_experimental);
N_unmeasured = size(XY,1) - num_new;

Sa = S(measured_indices_just_experimental,measured_indices_just_experimental);
Ra = C(measured_indices_just_experimental,:);
R = C;

jj = 1;
risk_ = []; cluster_path = []; query_path = []; certainty_path = {};

% while unmeasured samples still exist
fprintf('Number to measure: %d', N_unmeasured);
BO_lambda = .1;
cluster_similarity_matrix = NaN( size(XY,1), size(XY,1));
query_in_sequence = coarse_indices(:)';
query_BO_bool = false;
query_UCB_list = [];
query_AL_list = [];
for ii = 1:length(coarse_indices)
    query_list(ii).value = coarse_indices(ii);
    query_list(ii).type = 'initial';
end
AL_number_of_iterations = 0;
BO_number_of_iterations = 1;
fmiw_ = [];
Nl = 0;

% Active learning sequential sampling
while N_unmeasured > 0
    experiment.measured_indices_just_experimental = measured_indices_just_experimental;
    AL_number_of_iterations = AL_number_of_iterations + 1;
    brk_str = '\b\b\b';
    brk_length = length( num2str( N_unmeasured + 1) );
    fprintf([brk_str(1:brk_length*2) '%d'], N_unmeasured);
    
    % --------------- CLUSTERING --------------------------
    % Learn phase diagram
    [Ua, ~, ~, cluster_labels] = phase_diagram_model_AL_160718a(Xa, Ra, Sa, T, experiment);
    
    % ---------- PHASE REGION LABEL PROPAGATION ---------------
    [query, certainty, risk, full_cluster_labels] = active_learning_for_phase_diagram_160718b(S, XY, R, measured_indices_just_experimental, Ua, cluster_labels);
    query = query-Nl;
    full_cluster_labels = split_cluster_graph(S, full_cluster_labels, R);
    
    cluster_path = [cluster_path; full_cluster_labels'];

    
    % ------------ IDENTIFY IF AL TO BO SWITCH SHOULD HAPPEN -------------
    % if we eventually want BO: Compute cluster similarity to past (and store).
    clusters_to_compare = max(1,AL_number_of_iterations-4):AL_number_of_iterations;
    fmiw = zeros(1,4);
    if AL_number_of_iterations >=5 && experiment.bayesian_optimization
        [clust_perf_mp] = clustering_performance_measures_170930a(full_cluster_labels((Nl+1):end)', cluster_path(end-4:end,(Nl+1):end));
        fmiw = extractfield(clust_perf_mp,'fowlkes_mallows_index')';
        cluster_similarity_matrix( AL_number_of_iterations, clusters_to_compare) = fmiw;
        fmiw_ = [fmiw_; fmiw(end-1)];
        
        % --- Is phase mapping converging?
        if (fmiw(1) >= experiment.phase_mapping_threshold) && (AL_number_of_iterations > 5)
            query_BO_bool = true;
        end
    end    
    
    % ------------ SAMPLING -------------------
    % We enter this section with query for active learning, then decide
    % whether to overwrite it with in_sequence, random, or BO.
    
    % for sequential measurements
    unmeasured_sample_list_just_exp = setdiff(1:size(XY,1), measured_indices_just_experimental);% + Nl;
    
    if isempty(unmeasured_sample_list_just_exp)
        break;
    end
    
    BO_alpha_ = [];
    functional_property_estimate = NaN( size(Cexp_XY,1),1);
    functional_property_std = NaN( size(Cexp_XY,1),1);    
    points_of_interest_for_experimental = [];
    if query_BO_bool

        % ---------- DOWN SELECT FOR REGIONS OF INTEREST ------------------------------
        full_cluster_labels_experimental_points = full_cluster_labels( (Nl+1):end);
        unique_labels = unique(full_cluster_labels_experimental_points );
        max_for_each_region = 0*unique_labels;
        
        switch experiment.phase_region_condition
            % --------- COMPUTE GP FOR EACH REGION --------------
            case 'Phase_regions_with_max_by_GP'
                
                % Rather than calculating for all regions, can calculate just for
                % the regions of interest. But then regions of interest must be identified prior.
                [current_cluster_unmeasured_indices_list, current_cluster_measured_indices_list, current_cluster_all_indices_list, current_cluster_all_indices_func_prop_est_list, current_cluster_all_indices_func_prop_std_list] = deal(cell(length(unique_labels),1));
                current_cluster_label_list = cell(length(unique_labels),1);
                label_index = 1;
                for uu = unique_labels'
                    points_in_current_cluster = full_cluster_labels_experimental_points == uu;
                    current_cluster_measured_indices = measured_indices_just_experimental( full_cluster_labels_experimental_points(measured_indices_just_experimental) == uu );

                    current_cluster_unmeasured_indices_list{label_index} = unmeasured_sample_list_just_exp( full_cluster_labels_experimental_points(unmeasured_sample_list_just_exp) == uu );
                    current_cluster_measured_indices_list{label_index} = measured_indices_just_experimental( full_cluster_labels_experimental_points(measured_indices_just_experimental) == uu );
                    current_cluster_label_list{label_index} = uu;
                    current_cluster_all_indices_list{label_index} = find(points_in_current_cluster);
                    
                    if ~isempty(current_cluster_measured_indices)
                        gprMdl = fitrgp(Cexp_XY(current_cluster_measured_indices,:), ...
                            experiment.functional_property(current_cluster_measured_indices,:) );
                        [ypred,ystd] = predict(gprMdl,Cexp_XY(points_in_current_cluster,:));
                        functional_property_estimate(points_in_current_cluster) = ypred;
                        functional_property_std(points_in_current_cluster) = ystd;
                        
                        max_for_each_region(label_index) = max(ypred);
                        current_cluster_all_indices_func_prop_est_list{label_index} = ypred;
                        current_cluster_all_indices_func_prop_std_list{label_index} = ystd;
                    end
                    label_index = label_index + 1;
                end
                phase_region_table = table(current_cluster_label_list, max_for_each_region, current_cluster_all_indices_list, current_cluster_all_indices_func_prop_est_list, current_cluster_all_indices_func_prop_std_list, ...
                    current_cluster_measured_indices_list, current_cluster_unmeasured_indices_list);
                phase_region_table.Properties.VariableNames = {'label','max', 'all_idx', 'GPest', 'GPstd', 'measured_idx', 'unmeasured_idx'};
                phase_region_table = sortrows(phase_region_table, 'max', 'descend');
                
        end       
        
        % ---------- select region to investigate -------------
        % ------ this should take into account if more than 1 region to
        % search and no merging -----------
        region_to_investigate = cellfun(@isempty, phase_region_table.unmeasured_idx);
        region_to_investigate = find( ~region_to_investigate);
        region_to_investigate = region_to_investigate(1);
        all_points_of_interest = phase_region_table.all_idx{region_to_investigate};
        unmeasured_points_of_interest = phase_region_table.unmeasured_idx{region_to_investigate};
        
        % --------- COMPUTE UCB SAMPLING CRITERIA --------------
        Dsize = length(points_of_interest_for_experimental);
        
        % UCB with optional weight to focus on edges
        distance_to_edge = distance_to_edge_of_roi(Cexp_XY, all_points_of_interest, unmeasured_points_of_interest);
        BO_beta = 2*log(Dsize * (BO_number_of_iterations^2) * (pi^2) / (6* BO_lambda) );
        BO_alpha = functional_property_estimate(unmeasured_points_of_interest) ...
            + experiment.BO_exploration*sqrt(BO_beta) * functional_property_std(unmeasured_points_of_interest) ...
            - experiment.BO_edge_weight * distance_to_edge; % This is focus more on edges.
        
        % -------  point selection --------------------        
        [~, alpha_sorted_indices] = sort(BO_alpha, 'descend');
        sorted_unmeasured_points_of_interest = unmeasured_points_of_interest(alpha_sorted_indices);
                       
        query_BO = sorted_unmeasured_points_of_interest(1);
        query = query_BO;
        BO_number_of_iterations = BO_number_of_iterations + 1;
        BO_alpha_ = NaN(1,(size(C,1)-Nl));
        BO_alpha_(unmeasured_points_of_interest) = BO_alpha;
    end
    
    func_est = functional_property_estimate; func_std = functional_property_std;
    plot_active_learning_results(C, full_cluster_labels, certainty, fmiw_, query, query_list, query_BO_bool, func_est, func_std, BO_alpha_, Fa);

    
    % -------------- GETTING DATA FOR QUERY ------------------------
    % Function for measuring a new data point
    % query is only for experimental data (number doesn't include first indices for library entries).
    Xn = X(query,:);
    query_path = [query_path; query];
    
    if experiment.bayesian_optimization
        Fn = experiment.functional_property(query);
    end
    
    
    % -------------- DATA MANAGEMENT ------------------------------
    risk_ = [risk_ risk];
    certainty_path{jj} = certainty;
    F(jj) = jj; 
    jj = jj + 1;
    
    % add uncertain point
    if sum(measured_indices_just_experimental == query)
        '';
    end
    N_unmeasured = N_unmeasured - length(query);
    [measured_indices_just_experimental, sort_index] = sort([measured_indices_just_experimental query]);
    query_in_sequence = [query_in_sequence query];
    Fa_in_sequence = experiment.functional_property(query_in_sequence);
    Xa = [Xa; Xn];
    Xa = Xa([1:Nl Nl+sort_index], :);
    if experiment.bayesian_optimization
        Fa = [Fa; Fn];
        Fa = Fa(sort_index);
    end
    Ra = R([1:Nl Nl+measured_indices_just_experimental],:);
    Sa = S([1:Nl Nl+measured_indices_just_experimental], [1:Nl Nl+measured_indices_just_experimental]);
    
    if query_BO_bool
        query_UCB_list = [query_UCB_list query];
        query_list_length = length(query_list);
        query_list(query_list_length+1).value = query;
        query_list(query_list_length+1).type = 'UCB';
    else
        query_AL_list = [query_AL_list query];
        query_list_length = length(query_list);
        query_list(query_list_length+1).value = query;
        query_list(query_list_length+1).type = 'AL';        
    end
    
    if query_BO_bool && sum(max(experiment.functional_property) == Fa)
        break;
    end
    
    '.';
end

fprintf('\b0\n');

% results to pass
al_results.query_list = query_list;
al_results.cluster_path = cluster_path;
al_results.query_path = query_path;
al_results.certainty_path = certainty_path;
clearvars -except al_results;
'';

%------------- Use of Prior in S ----------------------
%------------------------------------------------------
function [S] = use_prior_to_impact_S(experiment, S)
if experiment.function_property_prior
    LL = 0;
    % W_boost is 1 if the prior labels are the same, 0 otherwise, except
    % for the diagonal which is 0.
    % Adding this will just boost samples with similar labels, but we also
    % want to subtract from samples with dissimilar labels.
    W_boost = experiment.W_boost;
    W_boost(W_boost == 0) = -1; % if no connection, set 0 -> -1
    W_boost = W_boost .* (~eye(size(W_boost,1))); % Connection =1, no connection = -1, diagonal = 0;
    W_boost_ = zeros(size(S));
    W_boost_((LL+1):end,(LL+1):end) = W_boost;
    S(S > 0) = 1 + experiment.S_boost_weight * W_boost_(S > 0);
end


%------------
function [distance_to_edge] = distance_to_edge_of_roi(XY_, all_points, unmeasured_points)

XY = XY_(all_points,:);
max_dist_from_center = max( pdist([mean(XY,1); XY]) );
distance_delta = .1;
if size(XY,1) > 2
    % create a voinoi diagram and identify the vertices for each point. 
    [V_poi,C_poi] = voronoin(XY);
    
    % select the vertices with a distance to the points_center that is over
    % 70% the max distance of the points from their mean.
    dist_vert_to_mean = squareform(pdist([mean(XY,1); V_poi]));
    dist_vert_to_mean = dist_vert_to_mean(2:end,1);
    vertex_idx = find( dist_vert_to_mean > .7 * max_dist_from_center );
    points_near_edge = false(size(XY,1),1);
    for i=1:length(C_poi)
        if ~isempty( intersect( C_poi{i}, vertex_idx  ) ); points_near_edge(i) = true; end
    end
    pne = find(points_near_edge);
    pnc = find(~points_near_edge);
    
    % points as the edge are given a distance of .1 from the edge.
    distance_to_edge = zeros( size(XY,1), 1 );
    distance_to_edge(pne) = distance_delta;
    
    % reorder the points: first edge points.
    XY_ordered = XY([ pne ; pnc ],: );
    
    % compute distance to the nearest edge point.
    D = squareform(pdist( XY_ordered ));
    D = D(1:length(pne),(length(pne)+1):end);
    min_dist_to_edge = min(D);
    distance_to_edge(pnc) = min_dist_to_edge + distance_delta;
    k = ismember(all_points, unmeasured_points);
    distance_to_edge = distance_to_edge(k);
else % only 1 point in set
    distance_to_edge = distance_delta * ones(length(unmeasured_points),1);
end

% ----------- Plot ----------------
% ----------------------------------
function plot_active_learning_results(C, full_cluster_labels, certainty, fmiw, query, query_path, BO_bool, func_est, func_std, BO_alpha, Fa)
figure(1); clf;
query_path = [query_path(:).value];
XY = tern2cart(C(:,[2 3 1]));
tx = [0 .6 .3 0]';
ty = [0 0 .3*sqrt(3) 0]';
x = XY(:,1); y = XY(:,2);

subplot(3,2,1);
terplot_general(tx, ty);
hold on;
scatter(x,y,(certainty*5).^2,full_cluster_labels,'filled');
plot(x(query),y(query),'rd','linewidth',2);
scatter(x(query_path),y(query_path),8^2,'k','linewidth',1);
hold off;

subplot(3,2,2);
plot(100*(1-fmiw));

if BO_bool
    subplot(3,2,3);
    terplot_general(tx, ty);
    hold on;
    yy = 1./(1 + exp(func_std) );
    scatter(x,y,yy*50,func_est,'filled');
    scatter(x(query_path),y(query_path),8^2,'k','linewidth',1);
    plot(x(query),y(query),'rd','linewidth',2);
    hold off;

    subplot(3,2,4);
    terplot_general(tx, ty);
    hold on;
    scatter(x,y,20,BO_alpha,'filled');
    scatter(x(query_path),y(query_path),8^2,'k','linewidth',1);
    plot(x(query),y(query),'rd','linewidth',2);
    hold off;

    BO_conv = (cummax(Fa(2:end)) - cummax(Fa(1:(end-1))));
    subplot(3,2,5);
    plot(BO_conv);
end
pause;

%------------ Load Data For Each Wafer -------------------
%---------------------------------------------------------
function [X, T, C, XY, clabels, experiment] = load_wafer_data(experiment)

R_type = 1;
X = []; T = []; C = []; XY = []; clabels = {}; library = [];
% if simulation, load appropriate simulation data

[X, T, C, XY, clabels, experiment] = load_FeGaPd(experiment);


experiment.R_type = R_type;

% ------- specific loading functions ---------------
%---------------------------------------------------
function [X, T, C, XY, clabels, experiment] = load_FeGaPd(experiment)
% experimental data
D1 = importdata('data/FeGaPd_CMP.txt');
D2 = importdata('data/FeGaPd_XRD.txt');
D3 = importdata('data/FeGaPd_Mag.txt');
labels_DFT = importdata('data/FeGaPd_DFT_regions.txt');

X = D2(2:end,:);
T = D2(1,:);
C = D1(1).data;
Mag = D3(:,1);
clabels = D1(1).colheaders;
XY = tern2cart(C(:,[2 3 1])*100);
% The W_boost does not use the modified Mag.
if experiment.phase_mapping_prior
    cluster_labels_func_props = labels_DFT;
    D_boost_data = squareform(pdist(cluster_labels_func_props));
    experiment.W_boost = ~D_boost_data - eye(size(X,1));    
end
if experiment.FeGaPd_Mag_boost_main_peak > 0
    max_points = Mag == 10;
    zz = (XY(:,1)-.19 ).^2 + (XY(:,2)- .05).^2;
    gg = 1*exp( -1*zz/.001);
    Boost_array = 0*Mag;
    Boost_array(max_points) = gg(max_points);
    Mag = Mag + Boost_array;
end

experiment.functional_property_modified = true;
experiment.functional_property = Mag;
