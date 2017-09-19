% Copyright (c) 2017 Zi Wang, Chengtao Li
% See also: gpopt.m
function results = batch_add_gpopt(objective, xmin, xmax, T, initx, inity, options)
% This function maximizes the function objective using BO with add-GP and 
% returns results as a cell of size 7, including the inferred argmax points 
% (guesses),the function values of the inferred argmax points (guessvals), the
% evaluated points (xx), the function values of the evaluated points
% (yy), the runtime to choose the points (choose_time) and extra time of
% inferring the argmax (extra_time).
% objective is a function handle;
% xmin is a column vector indicating the lower bound of the search space;
% xmax is a column vector indicating the upper bound of the search space;
% T is the number of sequential evaluations of the function;
% initx, inity are the initialization of observed inputs and their values.

if nargin <= 6
    options = struct();
end
if ~isfield(options, 'restart'); options.restart = 0; end
if ~isfield(options, 'bo_method'); options.bo_method = 'batch_ucb_dpp'; end
if ~isfield(options, 'batch_size'); options.batch_size = 10; end
if ~isfield(options, 'num_split'); options.num_split = 5; end
if ~isfield(options, 'savefilenm'); options.savefilenm = 'tmp.mat'; end
% When testing synthetic functions, one can add noise to the output.
if ~isfield(options, 'noiselevel'); options.noiselevel = 0; end
if ~isfield(options, 'seed'); options.seed = 42; end
% Set random seed
s = RandStream('mcg16807','Seed', options.seed);
RandStream.setGlobalStream(s);

% discretize each dimensions
options.discrete = zeros(length(xmin), options.num_split);
for i = 1:length(xmin)
    options.discrete(i,:) = linspace(xmin(i), xmax(i), options.num_split);
end

% Set options.restart = 1 to use the saved results and run more iterations
if options.restart && exist(savefilenm,'file') ~= 2
    options.restart = 0;
end
if options.restart == 0
    if isempty(initx)
        % initialize xx,yy with at least one pair of intx, inty
        initx = rand_sample_interval(xmin, xmax, 1);
        inity = objective(initx);
    end
    xx = initx;
    yy = inity;
    choose_time = []; % elapsed time to choose where to evaluate
    extra_time = []; % elapsed time to optimize mean function, hyper-parameters
    tstart = 0;
    
    [z, hyp] = sampleStructPriors(xx, yy, options);
    options = get_grid(z, options);
    options.z = z;

else
    restart_file = load(options.savefilenm);
    xx           =  restart_file.results{1};
    yy           =  restart_file.results{2};
    choose_time  =  restart_file.results{3};
    extra_time   =  restart_file.results{4};
    t            =  restart_file.results{5};
    z            =  restart_file.results{6};
    tstart = t;
    if tstart >= T
        return
    end
end

        
KernelMatrixInv = cell(1);

%% start optimization
for t = tstart+1 : T
    
    xnext = zeros(1, size(xx,2));
    
    if strcmp(options.bo_method, 'batch_rand')
        for batch_idx = 1:options.batch_size
            for i = 1:size(xx,2)
                xnext(i) = randsample(options.discrete(i,:),1);
            end
            xx = [ xx ; xnext ];
            yy = [ yy ; objective(xnext) + randn(1)*options.noiselevel];
        end
    else
        options = get_grid(z, options);
        
        tic
        
        % learn structure after every minibatch
        [z, hyp] = sampleStructPriors(xx, yy, options);
        options = get_grid(z, options);
        options.z = z;
        
        extra_time = [extra_time; toc];
        
        tic
        
        % Calculate and inverse the gram matrix
        KernelMatrix = compute_gram(xx, hyp, 1, z);
        KernelMatrixInv{1} = chol2invchol(KernelMatrix);

        all_cat = unique(z);

        % Start optimization group by group
        for i = 1:length(all_cat)
            coords = (z==all_cat(i));
            xx_sub = xx(:,coords);
            xmin_sub = xmin(coords);
            xmax_sub = xmax(coords);
            l = hyp.l(:,coords);
            sigma = hyp.sigma(:,all_cat(i));
            sigma0 = hyp.sigma0(:,all_cat(i));
            alpha = 1;
            beta = sqrt(size(xx_sub,2)*log(2*t)/5);
            optimum = ucb_choose(xx_sub, yy, KernelMatrixInv, [], ...
                sigma0, sigma, l, xmin_sub, xmax_sub, alpha, beta);

            xnext(coords) = optimum;
        end
        choose_time = [choose_time; toc];

        xx = [ xx ; xnext ];
        yy = [ yy ; objective(xnext) + rand(1) * options.noiselevel];
        
        % update inverse gram matrix
        KernelMatrix = compute_gram(xx, hyp, 1, z);
        KernelMatrixInv{1} = chol2invchol(KernelMatrix);
        
        xnextbatch = zeros(options.batch_size-1, size(xx,2));
        
        for i = 1:length(all_cat)
            coords = (z==all_cat(i));
            xx_sub = xx(:,coords);
            l = hyp.l(:,coords);
            sigma = hyp.sigma(:,all_cat(i));
            sigma0 = hyp.sigma0(:,all_cat(i));

            assert(length(find(z == all_cat(i))) <= 3);

            alpha = 1;
            beta = -sqrt(size(xx_sub,2)*log(2*t)/5);
            target_threshold = @(x) evaluateUCB(x, xx_sub, yy, KernelMatrixInv, l, sigma, sigma0, alpha, beta);
            beta = 2*sqrt(size(xx_sub,2)*log(2*t)/5);
            target = @(x) evaluateUCB(x, xx_sub, yy, KernelMatrixInv, l, sigma, sigma0, alpha, beta);

            options.curr_grid = filter_grid(options.xgrid{i}, target, target_threshold);            
            if size(options.curr_grid, 1) <= options.batch_size-1
                options.curr_grid = options.xgrid{i};
            end

            
            Kmm = computeKmm(options.curr_grid, l', sigma, sigma0);
            Kmn = computeKnm(options.curr_grid, xx_sub, l', sigma);
            Kmm = Kmm - Kmn * KernelMatrixInv{1} * Kmn';
            if size(Kmm, 1) <= options.batch_size-1
                C = 1:options.batch_size-1;
            elseif ~isempty(strfind(options.bo_method, 'batch_ucb_pe'))
                [C] = dpp_max(Kmm, options.batch_size-1);
            elseif ~isempty(strfind(options.bo_method, 'batch_ucb_dpp'))
                [C] = sample_dpp(decompose_kernel(Kmm), options.batch_size-1);
            else
                disp('Not Implemented!');
                pause;
            end
            xnextbatch(1:length(C),coords) = options.curr_grid(C,:);
        end
        
        if ~isempty(strfind(options.bo_method, 'fnc'))
            xnextbatchbyval = zeros(options.batch_size-1, size(xx,2));
            vals = zeros(options.batch_size-1, 1);

            for i = 1:length(all_cat)
                coords = (z==all_cat(i));
                xx_sub = xx(:,coords);
                l = hyp.l(:,coords);
                sigma = hyp.sigma(:,all_cat(i));
                sigma0 = hyp.sigma0(:,all_cat(i));
            
                assert(length(find(z == all_cat(i))) <= 3);

                alpha = 1;
                beta = sqrt(size(xx_sub,2)*log(2*t)/5);
                target = @(x) evaluateUCB(x, xx_sub, yy, KernelMatrixInv, l, sigma, sigma0, alpha, beta);

                curr_vals = zeros(options.batch_size-1, 1);
                for batch_idx = 1:options.batch_size-1
                    curr_vals(batch_idx) = target(xnextbatch(batch_idx,coords));
                end
                [~, new_idx] = sort(curr_vals, 'ascend');
                for batch_idx = 1:options.batch_size-1
                    xnextbatchbyval(batch_idx, coords) = xnextbatch(new_idx(batch_idx), coords);
                end
            end
            xnextbatch = xnextbatchbyval;
            
        elseif ~isempty(strfind(options.bo_method, 'ori'))
            xnextbatchbyval = zeros(options.batch_size-1, size(xx,2));
            vals = zeros(options.batch_size-1, 1);

            for i = 1:length(all_cat)
                coords = (z==all_cat(i));
                xx_sub = xx(:,coords);
                l = hyp.l(:,coords);
                sigma = hyp.sigma(:,all_cat(i));
                sigma0 = hyp.sigma0(:,all_cat(i));
                
                assert(length(find(z == all_cat(i))) <= 3);

                alpha = 1;
                beta = sqrt(size(xx_sub,2)*log(2*t)/5);
                target = @(x) evaluateUCB(x, xx_sub, yy, KernelMatrixInv, l, sigma, sigma0, alpha, beta);
                
                curr_vals = zeros(options.batch_size-1, 1);
                for batch_idx = 1:options.batch_size-1
                    curr_vals(batch_idx) = target(xnextbatch(batch_idx,coords));
                end
                
                curr_sums = repmat(curr_vals',options.batch_size-1,1) + repmat(vals,1,options.batch_size-1);
                [new_vals, new_idx] = sort(curr_sums(:), 'ascend');
                vals = new_vals(1:options.batch_size-1);
                new_idx = new_idx(1:options.batch_size-1);

                new_xnextbatchbyval = [];
                for batch_idx = 1:options.batch_size-1
                    tmp = xnextbatchbyval(mod((new_idx(batch_idx)-1), options.batch_size-1) + 1, :);
                    tmp(coords) = xnextbatch(floor((new_idx(batch_idx)-1) / (options.batch_size-1)) + 1, coords);
                    new_xnextbatchbyval = [new_xnextbatchbyval; tmp];
                end
                xnextbatchbyval = new_xnextbatchbyval; 
            end
            xnextbatch = xnextbatchbyval;
        end
        
        xx = [ xx; xnextbatch ];
        for i = 1:(options.batch_size-1)
            yy = [ yy; objective(xnextbatch(i,:)) + randn(1)*options.noiselevel];
        end
    end
        
    yy_to_show = max(yy);
    disp([num2str(t) ': ' 'val=' num2str(yy_to_show)])
    % save result every a few iterations
    if ~isempty(options.savefilenm) == 0
        results{1} = xx;
        results{2} = yy;
        results{3} = choose_time;
        results{4} = extra_time;
        results{5} = t;
        results{6} = z;
        options.savefilenm
        save(options.savefilenm, 'results');
    end
end
