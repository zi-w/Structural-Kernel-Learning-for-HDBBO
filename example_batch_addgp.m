% Copyright (c) 2017 Zi Wang, Chengtao Li
clear all; clc;
% add necessary paths
deploy_batch;

% Define function
dx = 20;
xmin = zeros(dx,1);
xmax = ones(dx,1);

f = sample_addGP(dx, dx, xmin, xmax);

% Save the file to a path 
options.savefilenm = 'tmp.mat';

% Set the GP hyper-parameters if you would like to fix them.
% Comment the following 3 lines out if you would like to learn them.
options.l = ones(1,dx)*50;
options.sigma = ones(1, dx)*5;
options.sigma0 = 0.0001*ones(1, dx);

% bo_method chosen from 
% 'batch_rand', {'batch_ucb_dpp', 'batch_ucb_pe'} * {'_ori', '_fnc'}
options.bo_method = 'batch_ucb_dpp_fnc';
% batch_size
options.batch_size = 10;
% split number for each dimension
options.num_split = 10;
% generate the grid of indices first
gen_grids(options.num_split);

% Start BO with batch add-GP
% The additive learning strategy is based on the paper
% Wang, Zi and Li, Chengtao and Jegelka, Stefanie and Kohli, Pushmeet. 
% Batched High-dimensional Bayesian Optimization via Structural Kernel 
% Learning. arXiv preprint arXiv:1703.01973.
batch_add_gpopt(f, xmin, xmax, 50, [], [], options)
