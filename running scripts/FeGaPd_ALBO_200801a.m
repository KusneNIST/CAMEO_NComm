function [] = FeGaPd_ALBO_200801a()

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
% Kusne, A. Gilad, et al. "On-the-fly closed-loop materials discovery
% via Bayesian active learning." Nature communications 11.1 (2020): 1-11.

% Follow the CAMEO algorithm.
% Parameters are set to be used with the Fe-Ga-Pd dataset.

% Path to needed folders and files
active_learning_directory = 'G:\My Drive\Research\AMRS\code\CAMEO submission 200817a\';

% save file name
tt = datestr(datetime,'yymmdd-HHMMss');
save_file_name = [active_learning_directory 'FGP_ALBO_' tt '.mat'];
cd(active_learning_directory);

% initialize
experiment.function_property_prior = false;
experiment.coarse_indices = randi(278);
experiment.graph_NN = 10;

% phase mapping
experiment.phase_mapping_prior = true;
experiment.num_endmembers = 3; % number of endmembers for GRENDEL
experiment.model_type = 1;
experiment.max_iter = 100; % limit GRENDEL to 100 iterations for each run.
experiment.GRENDEL_lower_bound_cutoff = .9;
experiment.graph_distance_multiplier = 1.2;   
experiment.distance_type = 'cosine';
experiment.sim_sigma = 1;
experiment.graphcut_balance = 0.75; 
experiment.graphcut_weight = 100; 
experiment.num_clusters = 5; % initialize number of clusters     
experiment.S_boost_weight = 0.5;

% Bayesian optimization parameters
experiment.BO_edge_weight = 10;
experiment.bayesian_optimization = true;
experiment.phase_region_condition = 'Phase_regions_with_max_by_GP';
experiment.BO_number_of_regions_to_search = 1;
experiment.BO_exploration = 1; % this multiples sqrt(beta)
experiment.phase_mapping_threshold = 0.8;
% Modify Remnant Magnetization for optimization.
experiment.FeGaPd_Mag_boost_main_peak = true; % main peak is flattened due to saturation, this smoothes the peak.

run_iteration.experiment = experiment;
run_iteration.results = run_active_learning_experiment_190618a2(experiment);

save(save_file_name);
