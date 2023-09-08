%% Read in experimental data
tStart = tic;
sheet_one = 'P3K1.xlsx';
sheet_two = 'P3K2.xlsx';
folder_name = 'P3K_7params_large_v3_testing_run ' + string(datetime('now'));
data1 = readmatrix(sheet_one);
data1 = data1.*0.0313;        % convert from AU to uM
% Repeat with replicate
data2 = readmatrix(sheet_two);
data2 = data2.*0.0313;        % convert from AU to uM
%%
%% Sample cells
rng('default');
cell_sample_size = 100;          % change this if changing cell sample size (KEEP IT AN EVEN NUMBER)
row_indices1 = randsample(1:length(data1), cell_sample_size/2);
rand_data1 = data1(row_indices1,:);
% Repeat with replicate
row_indices2 = randsample(1:length(data2), cell_sample_size/2);
rand_data2 = data2(row_indices2,:);
% Concatenate replicates
row_indices_all = [row_indices1'; row_indices2'];
rand_data_all = [rand_data1; rand_data2];
%% Sample more different cells
weights1 = ones(1,length(data1));
weights1(1,row_indices1) = 0;
weights2 = ones(1,length(data2));
weights2(1,row_indices2) = 0;
row_indices12 = datasample(1:length(data1), cell_sample_size/2, 'Weights', weights1, 'Replace', false);
row_indices22 = datasample(1:length(data2), cell_sample_size/2, 'Weights', weights2, 'Replace', false);
row_indices_all2 = [row_indices12'; row_indices22'];
rand_data12 = data1(row_indices12,:);
rand_data22 = data2(row_indices22,:);
rand_data_all2 = [rand_data12; rand_data22];
%% Sample more different cells
weights1 = ones(1,length(data1));
weights1(1,[row_indices1 row_indices12]) = 0;
weights2 = ones(1,length(data2));
weights2(1,[row_indices2 row_indices22]) = 0;
row_indices1 = datasample(1:length(data1), cell_sample_size/2, 'Weights', weights1, 'Replace', false);
row_indices2 = datasample(1:length(data2), cell_sample_size/2, 'Weights', weights2, 'Replace', false);
row_indices = [row_indices1'; row_indices2'];
rand_data1 = data1(row_indices1,:);
rand_data2 = data2(row_indices2,:);
rand_data = [rand_data1; rand_data2];
%%
%% Sample initial conditions
format short g      % displays more decimals
init_sample_size = 100;           % change this if changing number of starting initial conditions
ref = [0.1775, log(2)/9, 0.5188, 1E-6, 0.004, 14, 0.08];     % reference point (from 2021 paper)
x0 = [0.1775, log(2)/9, 0.5188, 1E-6, 0.004];       % parameters with 0.1x-10x constraint
delay0 = 14;
nfkb_init0 = 0.08;
% Get random number [-1,1] as power of 10 to then multiply reference point
rand_x = 10.^(-1 + (1-(-1)).*rand(init_sample_size,length(x0))).*x0;
% rand_x = log10(rand_x);
% Custom constraints
rand_delays = 7 + (28-7).*rand(init_sample_size,1);
rand_nfkb_inits = 0.04 + (0.3-0.04).*rand(init_sample_size,1);
rand_inits = [rand_x rand_delays rand_nfkb_inits];
%%
%% Set up simulations
lb = [0.01775, 0.1*log(2)/9, 0.05188, 1E-7, 0.0004, 7, 0.04];     % lower bound
ub = [1.775, 10*log(2)/9, 5.188, 1E-5, 0.04, 28, 0.3];        % upper bound
options.Display = 'off';
options.PlotFcns = [];
options.MaxIter = 100;

params_fit = zeros(cell_sample_size, numel(ref), init_sample_size);
%%
%% Run the parameter fitting algorithm
first_half_c = cell_sample_size/2;
sec_half_c = cell_sample_size/2+1;
mycluster = parcluster;
mycluster.NumWorkers = 50;
mycluster.parpool(50);
parfor i = 1:size(rand_inits,1)
   for c = 1:first_half_c
       xi = invert_constraint(rand_inits(i,:),lb,ub);
       fun = @(p) opt_rmsd(p, sheet_one, row_indices1(c)+1);
       [p_fit, fval, exitflag, output] = fminsearch(fun, xi, options);
       p_fit = constrain(p_fit,lb,ub);
       params_fit(c,:,i) = p_fit;
   end
   for c = sec_half_c:cell_sample_size
       xi = invert_constraint(rand_inits(i,:),lb,ub);
       fun = @(p) opt_rmsd(p, sheet_two, row_indices2(c-cell_sample_size/2)+1);
       [p_fit, fval, exitflag, output] = fminsearch(fun, xi, options);
       p_fit = constrain(p_fit,lb,ub);
       params_fit(c,:,i) = p_fit;
   end
   disp("parallel iteration " + string(i) + " completed.")
end
%%
%% Re-simulate with new parameters
end_rmsd = zeros(cell_sample_size, init_sample_size);
raw_rmsd = zeros(cell_sample_size, init_sample_size);
top_rmsds = zeros(cell_sample_size, 10);
best_params = zeros(cell_sample_size, numel(ref), 10);
mkdir(folder_name);
parfor c = 1:cell_sample_size
   sim_data = zeros(init_sample_size, 98);
   cur_end_rmsd = zeros(1,init_sample_size);
   for i = 1:init_sample_size
       p_mod = [6 3 params_fit(c,1,i); 15 1 params_fit(c,2,i); 16 1 params_fit(c,2,i); 66 1 params_fit(c,3,i); 68 1 params_fit(c,4,i); 75 1 params_fit(c,5,i); 6 4 params_fit(c,6,i)];
       init_mod = {'NFkB', params_fit(c,7,i)};
       doses = 100;
       dose_scale = 1/1500; % Convert to uM. Pam3CSK molecular weight: 1.5KDa
       names = {'NFkBn', 'IkBaNFkBn'};
       options = struct;
       options.DEBUG = 0;
       options.SIM_TIME = 8*60 + 5;
       % Simulate (only need to equilibrate on first iteration)
       output = [];
       if isempty(output)
           try
               [t,x,simdata2] = nfkbSimulate({'Pam3CSK',doses*dose_scale},names, p_mod, init_mod,options);
           catch
               warning('on');
               warning('nfkbSimulate failed. Proceeding to set rmsd to nan.');
               warning('off');
               cur_end_rmsd(i) = NaN;
               end_rmsd(c,i) = cur_end_rmsd(i);
               raw_rmsd(c,i) = cur_end_rmsd(i);
               continue
           end
       else
           options.STEADY_STATE = simdata2.STEADY_STATE;
           [~,x] = nfkbSimulate({'Pam3CSK',doses*dose_scale},names, p_mod, {},options);
       end
       output = cat(3,output,x);
       sim_nfkb_curve = squeeze(output(:,strcmp(names,'NFkBn'),:));
       sim_ikba_curve = squeeze(output(:,strcmp(names,'IkBaNFkBn'),:));
       sim_nfkb_curve = sim_nfkb_curve + sim_ikba_curve;
       sim_nfkb_curve = sim_nfkb_curve(1 : 5 : end);
       % bring sim to same zero point as experimental data
       zero = sim_nfkb_curve(1,1) - rand_data(c,1);
       sim_nfkb_curve = sim_nfkb_curve - zero;
       sim_data(i,:) = sim_nfkb_curve;
       dif = sim_nfkb_curve - rand_data(c,:)';
       raw_rmsd(c,i) = rms(dif, 'omitnan');
       dif(1:48,1) = dif(1:48,1).*sqrt(2);
       ratio = params_fit(c,:,i)./ref;
       ratio(1,1:5) = log10(ratio(1,1:5));
       ratio(1,6:7) = log2(ratio(1,6:7));
       cur_end_rmsd(i) = rms(dif, 'omitnan') + 0.1*rms(ratio);
       end_rmsd(c,i) = cur_end_rmsd(i);
   end
   % Get top 10% rmsds
   [cur_end_rmsd,indices] = sort(cur_end_rmsd);
   top_rmsds(c,:) = cur_end_rmsd(1:10);
   % end_rmsd(c,:) = end_rmsd(c,indices);
   top_sim_data = sim_data(indices',:);
   top_sim_data = top_sim_data(1:10,:);
   top_params = params_fit(c,:,indices);
   top_params = top_params(1,:,1:10);
   best_params(c,:,:) = top_params;
   % Plot data
   fig = figure;
   sums = sum(sim_data,2,'omitnan');
   [~,I] = sort(sums, 'descend');
   ranked_sim = sim_data(I',:);
   plot(top_sim_data')
   hold on
   % Ax = gca;
   % Ax.XDisplayLabels = nan(size(Ax.XDisplayData));
   % Ax.YDisplayLabels = nan(size(Ax.YDisplayData));
   % colormap default
   % Plot experimental data
   plot(rand_data(c,:)', 'LineWIdth', 3)
   if c <= first_half_c
       plot_name = 'P3Ki0' + string(1) + '_row_#' + string(row_indices(c)+1);
   else
       plot_name = 'P3Ki0' + string(2) + '_row_#' + string(row_indices(c)+1);
   end
   print(fig, plot_name, '-dpng')
   plot_file = plot_name + '.png';
   movefile(plot_file, folder_name);
   hold off
end
% Plot data
% sums = sum(rand_data,2,'omitnan');
% [~,I] = sort(sums, 'descend');
% ranked_exp = rand_data(I',:);
% plot(rand_data')
% Ax = gca;
% Ax.XDisplayLabels = nan(size(Ax.XDisplayData));
% Ax.YDisplayLabels = nan(size(Ax.YDisplayData));
% Ax.Title = 'Experimental';
% colormap default
%%
%% Store relevant data
result = cell(1,7);      % (1) cells, (2) initial conditions, (3) param fits
result{1} = [row_indices1 + 1; row_indices2 + 1];
result{2} = rand_inits;
result{3} = params_fit;
result{4} = raw_rmsd;
result{5} = end_rmsd;
result{6} = top_rmsds;
result{7} = best_params;
save('result_cell', 'result');
tEnd = toc(tStart);
save('runtime', 'tEnd');
movefile('result_cell.mat', folder_name);
movefile('runtime.mat', folder_name);
       
% NEW STATE
       
%% Read in experimental data
tStart = tic;
sheet_one = 'P3Kib1.xlsx';
sheet_two = 'P3Kib2.xlsx';
folder_name = 'P3Kib_7params_large_v3_testing_run ' + string(datetime('now'));
data1 = readmatrix(sheet_one);
data1 = data1.*0.0313;        % convert from AU to uM
% Repeat with replicate
data2 = readmatrix(sheet_two);
data2 = data2.*0.0313;        % convert from AU to uM
%%
%% Sample cells
rng('default');
cell_sample_size = 100;          % change this if changing cell sample size (KEEP IT AN EVEN NUMBER)
row_indices1 = randsample(1:length(data1), cell_sample_size/2);
rand_data1 = data1(row_indices1,:);
% Repeat with replicate
row_indices2 = randsample(1:length(data2), cell_sample_size/2);
rand_data2 = data2(row_indices2,:);
% Concatenate replicates
row_indices_all = [row_indices1'; row_indices2'];
rand_data_all = [rand_data1; rand_data2];
%% Sample more different cells
weights1 = ones(1,length(data1));
weights1(1,row_indices1) = 0;
weights2 = ones(1,length(data2));
weights2(1,row_indices2) = 0;
row_indices12 = datasample(1:length(data1), cell_sample_size/2, 'Weights', weights1, 'Replace', false);
row_indices22 = datasample(1:length(data2), cell_sample_size/2, 'Weights', weights2, 'Replace', false);
row_indices_all2 = [row_indices12'; row_indices22'];
rand_data12 = data1(row_indices12,:);
rand_data22 = data2(row_indices22,:);
rand_data_all2 = [rand_data12; rand_data22];
%% Sample more different cells
weights1 = ones(1,length(data1));
weights1(1,[row_indices1 row_indices12]) = 0;
weights2 = ones(1,length(data2));
weights2(1,[row_indices2 row_indices22]) = 0;
row_indices1 = datasample(1:length(data1), cell_sample_size/2, 'Weights', weights1, 'Replace', false);
row_indices2 = datasample(1:length(data2), cell_sample_size/2, 'Weights', weights2, 'Replace', false);
row_indices = [row_indices1'; row_indices2'];
rand_data1 = data1(row_indices1,:);
rand_data2 = data2(row_indices2,:);
rand_data = [rand_data1; rand_data2];
%%
%% Sample initial conditions
format short g      % displays more decimals
init_sample_size = 100;           % change this if changing number of starting initial conditions
ref = [0.1775, log(2)/9, 0.5188, 1E-6, 0.004, 14, 0.08];     % reference point (from 2021 paper)
x0 = [0.1775, log(2)/9, 0.5188, 1E-6, 0.004];       % parameters with 0.1x-10x constraint
delay0 = 14;
nfkb_init0 = 0.08;
% Get random number [-1,1] as power of 10 to then multiply reference point
rand_x = 10.^(-1 + (1-(-1)).*rand(init_sample_size,length(x0))).*x0;
% rand_x = log10(rand_x);
% Custom constraints
rand_delays = 7 + (28-7).*rand(init_sample_size,1);
rand_nfkb_inits = 0.04 + (0.3-0.04).*rand(init_sample_size,1);
rand_inits = [rand_x rand_delays rand_nfkb_inits];
%%
%% Set up simulations
lb = [0.01775, 0.1*log(2)/9, 0.05188, 1E-7, 0.0004, 7, 0.04];     % lower bound
ub = [1.775, 10*log(2)/9, 5.188, 1E-5, 0.04, 28, 0.3];        % upper bound
options.Display = 'off';
options.PlotFcns = [];
options.MaxIter = 100;

params_fit = zeros(cell_sample_size, numel(ref), init_sample_size);
%%
%% Run the parameter fitting algorithm
first_half_c = cell_sample_size/2;
sec_half_c = cell_sample_size/2+1;
% mycluster = parcluster;
% mycluster.NumWorkers = 50;
% mycluster.parpool(50);
parfor i = 1:size(rand_inits,1)
   for c = 1:first_half_c
       xi = invert_constraint(rand_inits(i,:),lb,ub);
       fun = @(p) opt_rmsd(p, sheet_one, row_indices1(c)+1);
       [p_fit, fval, exitflag, output] = fminsearch(fun, xi, options);
       p_fit = constrain(p_fit,lb,ub);
       params_fit(c,:,i) = p_fit;
   end
   for c = sec_half_c:cell_sample_size
       xi = invert_constraint(rand_inits(i,:),lb,ub);
       fun = @(p) opt_rmsd(p, sheet_two, row_indices2(c-cell_sample_size/2)+1);
       [p_fit, fval, exitflag, output] = fminsearch(fun, xi, options);
       p_fit = constrain(p_fit,lb,ub);
       params_fit(c,:,i) = p_fit;
   end
   disp("parallel iteration " + string(i) + " completed.")
end
%%
%% Re-simulate with new parameters
end_rmsd = zeros(cell_sample_size, init_sample_size);
raw_rmsd = zeros(cell_sample_size, init_sample_size);
top_rmsds = zeros(cell_sample_size, 10);
best_params = zeros(cell_sample_size, numel(ref), 10);
mkdir(folder_name);
parfor c = 1:cell_sample_size
   sim_data = zeros(init_sample_size, 98);
   cur_end_rmsd = zeros(1,init_sample_size);
   for i = 1:init_sample_size
       p_mod = [6 3 params_fit(c,1,i); 15 1 params_fit(c,2,i); 16 1 params_fit(c,2,i); 66 1 params_fit(c,3,i); 68 1 params_fit(c,4,i); 75 1 params_fit(c,5,i); 6 4 params_fit(c,6,i)];
       init_mod = {'NFkB', params_fit(c,7,i)};
       doses = 100;
       dose_scale = 1/1500; % Convert to uM. Pam3CSK molecular weight: 1.5KDa
       names = {'NFkBn', 'IkBaNFkBn'};
       options = struct;
       options.DEBUG = 0;
       options.SIM_TIME = 8*60 + 5;
       % Simulate (only need to equilibrate on first iteration)
       output = [];
       if isempty(output)
           try
               [t,x,simdata2] = nfkbSimulate({'Pam3CSK',doses*dose_scale},names, p_mod, init_mod,options);
           catch
               warning('on');
               warning('nfkbSimulate failed. Proceeding to set rmsd to nan.');
               warning('off');
               cur_end_rmsd(i) = NaN;
               end_rmsd(c,i) = cur_end_rmsd(i);
               raw_rmsd(c,i) = cur_end_rmsd(i);
               continue
           end
       else
           options.STEADY_STATE = simdata2.STEADY_STATE;
           [~,x] = nfkbSimulate({'Pam3CSK',doses*dose_scale},names, p_mod, {},options);
       end
       output = cat(3,output,x);
       sim_nfkb_curve = squeeze(output(:,strcmp(names,'NFkBn'),:));
       sim_ikba_curve = squeeze(output(:,strcmp(names,'IkBaNFkBn'),:));
       sim_nfkb_curve = sim_nfkb_curve + sim_ikba_curve;
       sim_nfkb_curve = sim_nfkb_curve(1 : 5 : end);
       % bring sim to same zero point as experimental data
       zero = sim_nfkb_curve(1,1) - rand_data(c,1);
       sim_nfkb_curve = sim_nfkb_curve - zero;
       sim_data(i,:) = sim_nfkb_curve;
       dif = sim_nfkb_curve - rand_data(c,:)';
       raw_rmsd(c,i) = rms(dif, 'omitnan');
       dif(1:48,1) = dif(1:48,1).*sqrt(2);
       ratio = params_fit(c,:,i)./ref;
       ratio(1,1:5) = log10(ratio(1,1:5));
       ratio(1,6:7) = log2(ratio(1,6:7));
       cur_end_rmsd(i) = rms(dif, 'omitnan') + 0.1*rms(ratio);
       end_rmsd(c,i) = cur_end_rmsd(i);
   end
   % Get top 10% rmsds
   [cur_end_rmsd,indices] = sort(cur_end_rmsd);
   top_rmsds(c,:) = cur_end_rmsd(1:10);
   top_sim_data = sim_data(indices',:);
   top_sim_data = top_sim_data(1:10,:);
   top_params = params_fit(c,:,indices);
   top_params = top_params(1,:,1:10);
   best_params(c,:,:) = top_params;
   % Plot data
   fig = figure;
   sums = sum(sim_data,2,'omitnan');
   [~,I] = sort(sums, 'descend');
   ranked_sim = sim_data(I',:);
   plot(top_sim_data')
   hold on
   % Ax = gca;
   % Ax.XDisplayLabels = nan(size(Ax.XDisplayData));
   % Ax.YDisplayLabels = nan(size(Ax.YDisplayData));
   % colormap default
   % Plot experimental data
   plot(rand_data(c,:)', 'LineWIdth', 3)
   if c <= first_half_c
       plot_name = 'P3Kib' + string(1) + '_row_#' + string(row_indices(c)+1);
   else
       plot_name = 'P3Kib' + string(2) + '_row_#' + string(row_indices(c)+1);
   end
   print(fig, plot_name, '-dpng')
   plot_file = plot_name + '.png';
   movefile(plot_file, folder_name);
   hold off
end
% Plot data
% sums = sum(rand_data,2,'omitnan');
% [~,I] = sort(sums, 'descend');
% ranked_exp = rand_data(I',:);
% plot(rand_data')
% Ax = gca;
% Ax.XDisplayLabels = nan(size(Ax.XDisplayData));
% Ax.YDisplayLabels = nan(size(Ax.YDisplayData));
% Ax.Title = 'Experimental';
% colormap default
%%
%% Store relevant data
result = cell(1,7);      % (1) cells, (2) initial conditions, (3) param fits
result{1} = [row_indices1 + 1; row_indices2 + 1];
result{2} = rand_inits;
result{3} = params_fit;
result{4} = raw_rmsd;
result{5} = end_rmsd;
result{6} = top_rmsds;
result{7} = best_params;
save('result_cell', 'result');
tEnd = toc(tStart);
save('runtime', 'tEnd');
movefile('result_cell.mat', folder_name);
movefile('runtime.mat', folder_name);
       
% NEW STATE
      
%% Read in experimental data
tStart = tic;
sheet_one = 'P3Kig1.xlsx';
sheet_two = 'P3Kig2.xlsx';
folder_name = 'P3Kig_7params_large_v3_testing_run ' + string(datetime('now'));
data1 = readmatrix(sheet_one);
data1 = data1.*0.0313;        % convert from AU to uM
% Repeat with replicate
data2 = readmatrix(sheet_two);
data2 = data2.*0.0313;        % convert from AU to uM
%%
%% Sample cells
rng('default');
cell_sample_size = 100;          % change this if changing cell sample size (KEEP IT AN EVEN NUMBER)
row_indices1 = randsample(1:length(data1), cell_sample_size/2);
rand_data1 = data1(row_indices1,:);
% Repeat with replicate
row_indices2 = randsample(1:length(data2), cell_sample_size/2);
rand_data2 = data2(row_indices2,:);
% Concatenate replicates
row_indices_all = [row_indices1'; row_indices2'];
rand_data_all = [rand_data1; rand_data2];
%% Sample more different cells
weights1 = ones(1,length(data1));
weights1(1,row_indices1) = 0;
weights2 = ones(1,length(data2));
weights2(1,row_indices2) = 0;
row_indices12 = datasample(1:length(data1), cell_sample_size/2, 'Weights', weights1, 'Replace', false);
row_indices22 = datasample(1:length(data2), cell_sample_size/2, 'Weights', weights2, 'Replace', false);
row_indices_all2 = [row_indices12'; row_indices22'];
rand_data12 = data1(row_indices12,:);
rand_data22 = data2(row_indices22,:);
rand_data_all2 = [rand_data12; rand_data22];
%% Sample more different cells
weights1 = ones(1,length(data1));
weights1(1,[row_indices1 row_indices12]) = 0;
weights2 = ones(1,length(data2));
weights2(1,[row_indices2 row_indices22]) = 0;
row_indices1 = datasample(1:length(data1), cell_sample_size/2, 'Weights', weights1, 'Replace', false);
row_indices2 = datasample(1:length(data2), cell_sample_size/2, 'Weights', weights2, 'Replace', false);
row_indices = [row_indices1'; row_indices2'];
rand_data1 = data1(row_indices1,:);
rand_data2 = data2(row_indices2,:);
rand_data = [rand_data1; rand_data2];
%%
%% Sample initial conditions
format short g      % displays more decimals
init_sample_size = 100;           % change this if changing number of starting initial conditions
ref = [0.1775, log(2)/9, 0.5188, 1E-6, 0.004, 14, 0.08];     % reference point (from 2021 paper)
x0 = [0.1775, log(2)/9, 0.5188, 1E-6, 0.004];       % parameters with 0.1x-10x constraint
delay0 = 14;
nfkb_init0 = 0.08;
% Get random number [-1,1] as power of 10 to then multiply reference point
rand_x = 10.^(-1 + (1-(-1)).*rand(init_sample_size,length(x0))).*x0;
% rand_x = log10(rand_x);
% Custom constraints
rand_delays = 7 + (28-7).*rand(init_sample_size,1);
rand_nfkb_inits = 0.04 + (0.3-0.04).*rand(init_sample_size,1);
rand_inits = [rand_x rand_delays rand_nfkb_inits];
%%
%% Set up simulations
lb = [0.01775, 0.1*log(2)/9, 0.05188, 1E-7, 0.0004, 7, 0.04];     % lower bound
ub = [1.775, 10*log(2)/9, 5.188, 1E-5, 0.04, 28, 0.3];        % upper bound
options.Display = 'off';
options.PlotFcns = [];
options.MaxIter = 100;

params_fit = zeros(cell_sample_size, numel(ref), init_sample_size);
%%
%% Run the parameter fitting algorithm
first_half_c = cell_sample_size/2;
sec_half_c = cell_sample_size/2+1;
% mycluster = parcluster;
% mycluster.NumWorkers = 50;
% mycluster.parpool(50);
parfor i = 1:size(rand_inits,1)
  for c = 1:first_half_c
      xi = invert_constraint(rand_inits(i,:),lb,ub);
      fun = @(p) opt_rmsd(p, sheet_one, row_indices1(c)+1);
      [p_fit, fval, exitflag, output] = fminsearch(fun, xi, options);
      p_fit = constrain(p_fit,lb,ub);
      params_fit(c,:,i) = p_fit;
  end
  for c = sec_half_c:cell_sample_size
      xi = invert_constraint(rand_inits(i,:),lb,ub);
      fun = @(p) opt_rmsd(p, sheet_two, row_indices2(c-cell_sample_size/2)+1);
      [p_fit, fval, exitflag, output] = fminsearch(fun, xi, options);
      p_fit = constrain(p_fit,lb,ub);
      params_fit(c,:,i) = p_fit;
  end
  disp("parallel iteration " + string(i) + " completed.")
end
%%
%% Re-simulate with new parameters
end_rmsd = zeros(cell_sample_size, init_sample_size);
raw_rmsd = zeros(cell_sample_size, init_sample_size);
top_rmsds = zeros(cell_sample_size, 10);
best_params = zeros(cell_sample_size, numel(ref), 10);
mkdir(folder_name);
parfor c = 1:cell_sample_size
  sim_data = zeros(init_sample_size, 98);
  cur_end_rmsd = zeros(1,init_sample_size);
  for i = 1:init_sample_size
      p_mod = [6 3 params_fit(c,1,i); 15 1 params_fit(c,2,i); 16 1 params_fit(c,2,i); 66 1 params_fit(c,3,i); 68 1 params_fit(c,4,i); 75 1 params_fit(c,5,i); 6 4 params_fit(c,6,i)];
      init_mod = {'NFkB', params_fit(c,7,i)};
      doses = 100;
      dose_scale = 1/1500; % Convert to uM. Pam3CSK molecular weight: 1.5KDa
      names = {'NFkBn', 'IkBaNFkBn'};
      options = struct;
      options.DEBUG = 0;
      options.SIM_TIME = 8*60 + 5;
      % Simulate (only need to equilibrate on first iteration)
      output = [];
      if isempty(output)
          try
              [t,x,simdata2] = nfkbSimulate({'Pam3CSK',doses*dose_scale},names, p_mod, init_mod,options);
          catch
              warning('on');
              warning('nfkbSimulate failed. Proceeding to set rmsd to nan.');
              warning('off');
              cur_end_rmsd(i) = NaN;
              end_rmsd(c,i) = cur_end_rmsd(i);
              raw_rmsd(c,i) = cur_end_rmsd(i);
              continue
          end
      else
          options.STEADY_STATE = simdata2.STEADY_STATE;
          [~,x] = nfkbSimulate({'Pam3CSK',doses*dose_scale},names, p_mod, {},options);
      end
      output = cat(3,output,x);
      sim_nfkb_curve = squeeze(output(:,strcmp(names,'NFkBn'),:));
      sim_ikba_curve = squeeze(output(:,strcmp(names,'IkBaNFkBn'),:));
      sim_nfkb_curve = sim_nfkb_curve + sim_ikba_curve;
      sim_nfkb_curve = sim_nfkb_curve(1 : 5 : end);
      % bring sim to same zero point as experimental data
      zero = sim_nfkb_curve(1,1) - rand_data(c,1);
      sim_nfkb_curve = sim_nfkb_curve - zero;
      sim_data(i,:) = sim_nfkb_curve;
      dif = sim_nfkb_curve - rand_data(c,:)';
      raw_rmsd(c,i) = rms(dif, 'omitnan');
      dif(1:48,1) = dif(1:48,1).*sqrt(2);
      ratio = params_fit(c,:,i)./ref;
      ratio(1,1:5) = log10(ratio(1,1:5));
      ratio(1,6:7) = log2(ratio(1,6:7));
      cur_end_rmsd(i) = rms(dif, 'omitnan') + 0.1*rms(ratio);
      end_rmsd(c,i) = cur_end_rmsd(i);
  end
  % Get top 10% rmsds
  [cur_end_rmsd,indices] = sort(cur_end_rmsd);
  top_rmsds(c,:) = cur_end_rmsd(1:10);
  top_sim_data = sim_data(indices',:);
  top_sim_data = top_sim_data(1:10,:);
  top_params = params_fit(c,:,indices);
  top_params = top_params(1,:,1:10);
  best_params(c,:,:) = top_params;
  % Plot data
  fig = figure;
  sums = sum(sim_data,2,'omitnan');
  [~,I] = sort(sums, 'descend');
  ranked_sim = sim_data(I',:);
  plot(top_sim_data')
  hold on
  % Ax = gca;
  % Ax.XDisplayLabels = nan(size(Ax.XDisplayData));
  % Ax.YDisplayLabels = nan(size(Ax.YDisplayData));
  % colormap default
  % Plot experimental data
  plot(rand_data(c,:)', 'LineWIdth', 3)
  if c <= first_half_c
      plot_name = 'P3Kig' + string(1) + '_row_#' + string(row_indices(c)+1);
  else
      plot_name = 'P3Kig' + string(2) + '_row_#' + string(row_indices(c)+1);
  end
  print(fig, plot_name, '-dpng')
  plot_file = plot_name + '.png';
  movefile(plot_file, folder_name);
  hold off
end
% Plot data
% sums = sum(rand_data,2,'omitnan');
% [~,I] = sort(sums, 'descend');
% ranked_exp = rand_data(I',:);
% plot(rand_data')
% Ax = gca;
% Ax.XDisplayLabels = nan(size(Ax.XDisplayData));
% Ax.YDisplayLabels = nan(size(Ax.YDisplayData));
% Ax.Title = 'Experimental';
% colormap default
%%
%% Store relevant data
result = cell(1,7);      % (1) cells, (2) initial conditions, (3) param fits
result{1} = [row_indices1 + 1; row_indices2 + 1];
result{2} = rand_inits;
result{3} = params_fit;
result{4} = raw_rmsd;
result{5} = end_rmsd;
result{6} = top_rmsds;
result{7} = best_params;
save('result_cell', 'result');
tEnd = toc(tStart);
save('runtime', 'tEnd');
movefile('result_cell.mat', folder_name);
movefile('runtime.mat', folder_name);

% NEW STATE
      
%% Read in experimental data
tStart = tic;
sheet_one = 'P3Ki31.xlsx';
sheet_two = 'P3Ki32.xlsx';
folder_name = 'P3Ki3_7params_large_v3_testing_run ' + string(datetime('now'));
data1 = readmatrix(sheet_one);
data1 = data1.*0.0313;        % convert from AU to uM
% Repeat with replicate
data2 = readmatrix(sheet_two);
data2 = data2.*0.0313;        % convert from AU to uM
%%
%% Sample cells
rng('default');
cell_sample_size = 100;          % change this if changing cell sample size (KEEP IT AN EVEN NUMBER)
row_indices1 = randsample(1:length(data1), cell_sample_size/2);
rand_data1 = data1(row_indices1,:);
% Repeat with replicate
row_indices2 = randsample(1:length(data2), cell_sample_size/2);
rand_data2 = data2(row_indices2,:);
% Concatenate replicates
row_indices_all = [row_indices1'; row_indices2'];
rand_data_all = [rand_data1; rand_data2];
%% Sample more different cells
weights1 = ones(1,length(data1));
weights1(1,row_indices1) = 0;
weights2 = ones(1,length(data2));
weights2(1,row_indices2) = 0;
row_indices12 = datasample(1:length(data1), cell_sample_size/2, 'Weights', weights1, 'Replace', false);
row_indices22 = datasample(1:length(data2), cell_sample_size/2, 'Weights', weights2, 'Replace', false);
row_indices_all2 = [row_indices12'; row_indices22'];
rand_data12 = data1(row_indices12,:);
rand_data22 = data2(row_indices22,:);
rand_data_all2 = [rand_data12; rand_data22];
%% Sample more different cells
weights1 = ones(1,length(data1));
weights1(1,[row_indices1 row_indices12]) = 0;
weights2 = ones(1,length(data2));
weights2(1,[row_indices2 row_indices22]) = 0;
row_indices1 = datasample(1:length(data1), cell_sample_size/2, 'Weights', weights1, 'Replace', false);
row_indices2 = datasample(1:length(data2), cell_sample_size/2, 'Weights', weights2, 'Replace', false);
row_indices = [row_indices1'; row_indices2'];
rand_data1 = data1(row_indices1,:);
rand_data2 = data2(row_indices2,:);
rand_data = [rand_data1; rand_data2];
%%
%% Sample initial conditions
format short g      % displays more decimals
init_sample_size = 100;           % change this if changing number of starting initial conditions
ref = [0.1775, log(2)/9, 0.5188, 1E-6, 0.004, 14, 0.08];     % reference point (from 2021 paper)
x0 = [0.1775, log(2)/9, 0.5188, 1E-6, 0.004];       % parameters with 0.1x-10x constraint
delay0 = 14;
nfkb_init0 = 0.08;
% Get random number [-1,1] as power of 10 to then multiply reference point
rand_x = 10.^(-1 + (1-(-1)).*rand(init_sample_size,length(x0))).*x0;
% rand_x = log10(rand_x);
% Custom constraints
rand_delays = 7 + (28-7).*rand(init_sample_size,1);
rand_nfkb_inits = 0.04 + (0.3-0.04).*rand(init_sample_size,1);
rand_inits = [rand_x rand_delays rand_nfkb_inits];
%%
%% Set up simulations
lb = [0.01775, 0.1*log(2)/9, 0.05188, 1E-7, 0.0004, 7, 0.04];     % lower bound
ub = [1.775, 10*log(2)/9, 5.188, 1E-5, 0.04, 28, 0.3];        % upper bound
options.Display = 'off';
options.PlotFcns = [];
options.MaxIter = 100;

params_fit = zeros(cell_sample_size, numel(ref), init_sample_size);
%%
%% Run the parameter fitting algorithm
first_half_c = cell_sample_size/2;
sec_half_c = cell_sample_size/2+1;
% mycluster = parcluster;
% mycluster.NumWorkers = 50;
% mycluster.parpool(50);
parfor i = 1:size(rand_inits,1)
  for c = 1:first_half_c
      xi = invert_constraint(rand_inits(i,:),lb,ub);
      fun = @(p) opt_rmsd(p, sheet_one, row_indices1(c)+1);
      [p_fit, fval, exitflag, output] = fminsearch(fun, xi, options);
      p_fit = constrain(p_fit,lb,ub);
      params_fit(c,:,i) = p_fit;
  end
  for c = sec_half_c:cell_sample_size
      xi = invert_constraint(rand_inits(i,:),lb,ub);
      fun = @(p) opt_rmsd(p, sheet_two, row_indices2(c-cell_sample_size/2)+1);
      [p_fit, fval, exitflag, output] = fminsearch(fun, xi, options);
      p_fit = constrain(p_fit,lb,ub);
      params_fit(c,:,i) = p_fit;
  end
  disp("parallel iteration " + string(i) + " completed.")
end
%%
%% Re-simulate with new parameters
end_rmsd = zeros(cell_sample_size, init_sample_size);
raw_rmsd = zeros(cell_sample_size, init_sample_size);
top_rmsds = zeros(cell_sample_size, 10);
best_params = zeros(cell_sample_size, numel(ref), 10);
mkdir(folder_name);
parfor c = 1:cell_sample_size
  sim_data = zeros(init_sample_size, 98);
  cur_end_rmsd = zeros(1,init_sample_size);
  for i = 1:init_sample_size
      p_mod = [6 3 params_fit(c,1,i); 15 1 params_fit(c,2,i); 16 1 params_fit(c,2,i); 66 1 params_fit(c,3,i); 68 1 params_fit(c,4,i); 75 1 params_fit(c,5,i); 6 4 params_fit(c,6,i)];
      init_mod = {'NFkB', params_fit(c,7,i)};
      doses = 100;
      dose_scale = 1/1500; % Convert to uM. Pam3CSK molecular weight: 1.5KDa
      names = {'NFkBn', 'IkBaNFkBn'};
      options = struct;
      options.DEBUG = 0;
      options.SIM_TIME = 8*60 + 5;
      % Simulate (only need to equilibrate on first iteration)
      output = [];
      if isempty(output)
          try
              [t,x,simdata2] = nfkbSimulate({'Pam3CSK',doses*dose_scale},names, p_mod, init_mod,options);
          catch
              warning('on');
              warning('nfkbSimulate failed. Proceeding to set rmsd to nan.');
              warning('off');
              cur_end_rmsd(i) = NaN;
              end_rmsd(c,i) = cur_end_rmsd(i);
              raw_rmsd(c,i) = cur_end_rmsd(i);
              continue
          end
      else
          options.STEADY_STATE = simdata2.STEADY_STATE;
          [~,x] = nfkbSimulate({'Pam3CSK',doses*dose_scale},names, p_mod, {},options);
      end
      output = cat(3,output,x);
      sim_nfkb_curve = squeeze(output(:,strcmp(names,'NFkBn'),:));
      sim_ikba_curve = squeeze(output(:,strcmp(names,'IkBaNFkBn'),:));
      sim_nfkb_curve = sim_nfkb_curve + sim_ikba_curve;
      sim_nfkb_curve = sim_nfkb_curve(1 : 5 : end);
      % bring sim to same zero point as experimental data
      zero = sim_nfkb_curve(1,1) - rand_data(c,1);
      sim_nfkb_curve = sim_nfkb_curve - zero;
      sim_data(i,:) = sim_nfkb_curve;
      dif = sim_nfkb_curve - rand_data(c,:)';
      raw_rmsd(c,i) = rms(dif, 'omitnan');
      dif(1:48,1) = dif(1:48,1).*sqrt(2);
      ratio = params_fit(c,:,i)./ref;
      ratio(1,1:5) = log10(ratio(1,1:5));
      ratio(1,6:7) = log2(ratio(1,6:7));
      cur_end_rmsd(i) = rms(dif, 'omitnan') + 0.1*rms(ratio);
      end_rmsd(c,i) = cur_end_rmsd(i);
  end
  % Get top 10% rmsds
  [cur_end_rmsd,indices] = sort(cur_end_rmsd);
  top_rmsds(c,:) = cur_end_rmsd(1:10);
  top_sim_data = sim_data(indices',:);
  top_sim_data = top_sim_data(1:10,:);
  top_params = params_fit(c,:,indices);
  top_params = top_params(1,:,1:10);
  best_params(c,:,:) = top_params;
  % Plot data
  fig = figure;
  sums = sum(sim_data,2,'omitnan');
  [~,I] = sort(sums, 'descend');
  ranked_sim = sim_data(I',:);
  plot(top_sim_data')
  hold on
  % Ax = gca;
  % Ax.XDisplayLabels = nan(size(Ax.XDisplayData));
  % Ax.YDisplayLabels = nan(size(Ax.YDisplayData));
  % colormap default
  % Plot experimental data
  plot(rand_data(c,:)', 'LineWIdth', 3)
  if c <= first_half_c
      plot_name = 'P3Ki3' + string(1) + '_row_#' + string(row_indices(c)+1);
  else
      plot_name = 'P3Ki3' + string(2) + '_row_#' + string(row_indices(c)+1);
  end
  print(fig, plot_name, '-dpng')
  plot_file = plot_name + '.png';
  movefile(plot_file, folder_name);
  hold off
end
% Plot data
% sums = sum(rand_data,2,'omitnan');
% [~,I] = sort(sums, 'descend');
% ranked_exp = rand_data(I',:);
% plot(rand_data')
% Ax = gca;
% Ax.XDisplayLabels = nan(size(Ax.XDisplayData));
% Ax.YDisplayLabels = nan(size(Ax.YDisplayData));
% Ax.Title = 'Experimental';
% colormap default
%%
%% Store relevant data
result = cell(1,7);      % (1) cells, (2) initial conditions, (3) param fits
result{1} = [row_indices1 + 1; row_indices2 + 1];
result{2} = rand_inits;
result{3} = params_fit;
result{4} = raw_rmsd;
result{5} = end_rmsd;
result{6} = top_rmsds;
result{7} = best_params;
save('result_cell', 'result');
tEnd = toc(tStart);
save('runtime', 'tEnd');
movefile('result_cell.mat', folder_name);
movefile('runtime.mat', folder_name);

% NEW STATE
      
%% Read in experimental data
tStart = tic;
sheet_one = 'P3Ki41.xlsx';
sheet_two = 'P3Ki42.xlsx';
folder_name = 'P3Ki4_7params_large_v3_testing_run ' + string(datetime('now'));
data1 = readmatrix(sheet_one);
data1 = data1.*0.0313;        % convert from AU to uM
% Repeat with replicate
data2 = readmatrix(sheet_two);
data2 = data2.*0.0313;        % convert from AU to uM
%%
%% Sample cells
rng('default');
cell_sample_size = 100;          % change this if changing cell sample size (KEEP IT AN EVEN NUMBER)
row_indices1 = randsample(1:length(data1), cell_sample_size/2);
rand_data1 = data1(row_indices1,:);
% Repeat with replicate
row_indices2 = randsample(1:length(data2), cell_sample_size/2);
rand_data2 = data2(row_indices2,:);
% Concatenate replicates
row_indices_all = [row_indices1'; row_indices2'];
rand_data_all = [rand_data1; rand_data2];
%% Sample more different cells
weights1 = ones(1,length(data1));
weights1(1,row_indices1) = 0;
weights2 = ones(1,length(data2));
weights2(1,row_indices2) = 0;
row_indices12 = datasample(1:length(data1), cell_sample_size/2, 'Weights', weights1, 'Replace', false);
row_indices22 = datasample(1:length(data2), cell_sample_size/2, 'Weights', weights2, 'Replace', false);
row_indices_all2 = [row_indices12'; row_indices22'];
rand_data12 = data1(row_indices12,:);
rand_data22 = data2(row_indices22,:);
rand_data_all2 = [rand_data12; rand_data22];
%% Sample more different cells
weights1 = ones(1,length(data1));
weights1(1,[row_indices1 row_indices12]) = 0;
weights2 = ones(1,length(data2));
weights2(1,[row_indices2 row_indices22]) = 0;
row_indices1 = datasample(1:length(data1), cell_sample_size/2, 'Weights', weights1, 'Replace', false);
row_indices2 = datasample(1:length(data2), cell_sample_size/2, 'Weights', weights2, 'Replace', false);
row_indices = [row_indices1'; row_indices2'];
rand_data1 = data1(row_indices1,:);
rand_data2 = data2(row_indices2,:);
rand_data = [rand_data1; rand_data2];
%%
%% Sample initial conditions
format short g      % displays more decimals
init_sample_size = 100;           % change this if changing number of starting initial conditions
ref = [0.1775, log(2)/9, 0.5188, 1E-6, 0.004, 14, 0.08];     % reference point (from 2021 paper)
x0 = [0.1775, log(2)/9, 0.5188, 1E-6, 0.004];       % parameters with 0.1x-10x constraint
delay0 = 14;
nfkb_init0 = 0.08;
% Get random number [-1,1] as power of 10 to then multiply reference point
rand_x = 10.^(-1 + (1-(-1)).*rand(init_sample_size,length(x0))).*x0;
% rand_x = log10(rand_x);
% Custom constraints
rand_delays = 7 + (28-7).*rand(init_sample_size,1);
rand_nfkb_inits = 0.04 + (0.3-0.04).*rand(init_sample_size,1);
rand_inits = [rand_x rand_delays rand_nfkb_inits];
%%
%% Set up simulations
lb = [0.01775, 0.1*log(2)/9, 0.05188, 1E-7, 0.0004, 7, 0.04];     % lower bound
ub = [1.775, 10*log(2)/9, 5.188, 1E-5, 0.04, 28, 0.3];        % upper bound
options.Display = 'off';
options.PlotFcns = [];
options.MaxIter = 100;

params_fit = zeros(cell_sample_size, numel(ref), init_sample_size);
%%
%% Run the parameter fitting algorithm
first_half_c = cell_sample_size/2;
sec_half_c = cell_sample_size/2+1;
% mycluster = parcluster;
% mycluster.NumWorkers = 50;
% mycluster.parpool(50);
parfor i = 1:size(rand_inits,1)
  for c = 1:first_half_c
      xi = invert_constraint(rand_inits(i,:),lb,ub);
      fun = @(p) opt_rmsd(p, sheet_one, row_indices1(c)+1);
      [p_fit, fval, exitflag, output] = fminsearch(fun, xi, options);
      p_fit = constrain(p_fit,lb,ub);
      params_fit(c,:,i) = p_fit;
  end
  for c = sec_half_c:cell_sample_size
      xi = invert_constraint(rand_inits(i,:),lb,ub);
      fun = @(p) opt_rmsd(p, sheet_two, row_indices2(c-cell_sample_size/2)+1);
      [p_fit, fval, exitflag, output] = fminsearch(fun, xi, options);
      p_fit = constrain(p_fit,lb,ub);
      params_fit(c,:,i) = p_fit;
  end
  disp("parallel iteration " + string(i) + " completed.")
end
%%
%% Re-simulate with new parameters
end_rmsd = zeros(cell_sample_size, init_sample_size);
raw_rmsd = zeros(cell_sample_size, init_sample_size);
top_rmsds = zeros(cell_sample_size, 10);
best_params = zeros(cell_sample_size, numel(ref), 10);
mkdir(folder_name);
parfor c = 1:cell_sample_size
  sim_data = zeros(init_sample_size, 98);
  cur_end_rmsd = zeros(1,init_sample_size);
  for i = 1:init_sample_size
      p_mod = [6 3 params_fit(c,1,i); 15 1 params_fit(c,2,i); 16 1 params_fit(c,2,i); 66 1 params_fit(c,3,i); 68 1 params_fit(c,4,i); 75 1 params_fit(c,5,i); 6 4 params_fit(c,6,i)];
      init_mod = {'NFkB', params_fit(c,7,i)};
      doses = 100;
      dose_scale = 1/1500; % Convert to uM. Pam3CSK molecular weight: 1.5KDa
      names = {'NFkBn', 'IkBaNFkBn'};
      options = struct;
      options.DEBUG = 0;
      options.SIM_TIME = 8*60 + 5;
      % Simulate (only need to equilibrate on first iteration)
      output = [];
      if isempty(output)
          try
              [t,x,simdata2] = nfkbSimulate({'Pam3CSK',doses*dose_scale},names, p_mod, init_mod,options);
          catch
              warning('on');
              warning('nfkbSimulate failed. Proceeding to set rmsd to nan.');
              warning('off');
              cur_end_rmsd(i) = NaN;
              end_rmsd(c,i) = cur_end_rmsd(i);
              raw_rmsd(c,i) = cur_end_rmsd(i);
              continue
          end
      else
          options.STEADY_STATE = simdata2.STEADY_STATE;
          [~,x] = nfkbSimulate({'Pam3CSK',doses*dose_scale},names, p_mod, {},options);
      end
      output = cat(3,output,x);
      sim_nfkb_curve = squeeze(output(:,strcmp(names,'NFkBn'),:));
      sim_ikba_curve = squeeze(output(:,strcmp(names,'IkBaNFkBn'),:));
      sim_nfkb_curve = sim_nfkb_curve + sim_ikba_curve;
      sim_nfkb_curve = sim_nfkb_curve(1 : 5 : end);
      % bring sim to same zero point as experimental data
      zero = sim_nfkb_curve(1,1) - rand_data(c,1);
      sim_nfkb_curve = sim_nfkb_curve - zero;
      sim_data(i,:) = sim_nfkb_curve;
      dif = sim_nfkb_curve - rand_data(c,:)';
      raw_rmsd(c,i) = rms(dif, 'omitnan');
      dif(1:48,1) = dif(1:48,1).*sqrt(2);
      ratio = params_fit(c,:,i)./ref;
      ratio(1,1:5) = log10(ratio(1,1:5));
      ratio(1,6:7) = log2(ratio(1,6:7));
      cur_end_rmsd(i) = rms(dif, 'omitnan') + 0.1*rms(ratio);
      end_rmsd(c,i) = cur_end_rmsd(i);
  end
  % Get top 10% rmsds
  [cur_end_rmsd,indices] = sort(cur_end_rmsd);
  top_rmsds(c,:) = cur_end_rmsd(1:10);
  top_sim_data = sim_data(indices',:);
  top_sim_data = top_sim_data(1:10,:);
  top_params = params_fit(c,:,indices);
  top_params = top_params(1,:,1:10);
  best_params(c,:,:) = top_params;
  % Plot data
  fig = figure;
  sums = sum(sim_data,2,'omitnan');
  [~,I] = sort(sums, 'descend');
  ranked_sim = sim_data(I',:);
  plot(top_sim_data')
  hold on
  % Ax = gca;
  % Ax.XDisplayLabels = nan(size(Ax.XDisplayData));
  % Ax.YDisplayLabels = nan(size(Ax.YDisplayData));
  % colormap default
  % Plot experimental data
  plot(rand_data(c,:)', 'LineWIdth', 3)
  if c <= first_half_c
      plot_name = 'P3Ki4' + string(1) + '_row_#' + string(row_indices(c)+1);
  else
      plot_name = 'P3Ki4' + string(2) + '_row_#' + string(row_indices(c)+1);
  end
  print(fig, plot_name, '-dpng')
  plot_file = plot_name + '.png';
  movefile(plot_file, folder_name);
  hold off
end
% Plot data
% sums = sum(rand_data,2,'omitnan');
% [~,I] = sort(sums, 'descend');
% ranked_exp = rand_data(I',:);
% plot(rand_data')
% Ax = gca;
% Ax.XDisplayLabels = nan(size(Ax.XDisplayData));
% Ax.YDisplayLabels = nan(size(Ax.YDisplayData));
% Ax.Title = 'Experimental';
% colormap default
%%
%% Store relevant data
result = cell(1,7);      % (1) cells, (2) initial conditions, (3) param fits
result{1} = [row_indices1 + 1; row_indices2 + 1];
result{2} = rand_inits;
result{3} = params_fit;
result{4} = raw_rmsd;
result{5} = end_rmsd;
result{6} = top_rmsds;
result{7} = best_params;
save('result_cell', 'result');
tEnd = toc(tStart);
save('runtime', 'tEnd');
movefile('result_cell.mat', folder_name);
movefile('runtime.mat', folder_name);
       
% NEW STATE
     
%% Read in experimental data
tStart = tic;
sheet_one = 'P3Ki01.xlsx';
sheet_two = 'P3Ki02.xlsx';
folder_name = 'P3Ki0_7params_large_v3_testing_run ' + string(datetime('now'));
data1 = readmatrix(sheet_one);
data1 = data1.*0.0313;        % convert from AU to uM
% Repeat with replicate
data2 = readmatrix(sheet_two);
data2 = data2.*0.0313;        % convert from AU to uM
%%
%% Sample cells
rng('default');
cell_sample_size = 100;          % change this if changing cell sample size (KEEP IT AN EVEN NUMBER)
row_indices1 = randsample(1:length(data1), cell_sample_size/2);
rand_data1 = data1(row_indices1,:);
% Repeat with replicate
row_indices2 = randsample(1:length(data2), cell_sample_size/2);
rand_data2 = data2(row_indices2,:);
% Concatenate replicates
row_indices_all = [row_indices1'; row_indices2'];
rand_data_all = [rand_data1; rand_data2];
%% Sample more different cells
weights1 = ones(1,length(data1));
weights1(1,row_indices1) = 0;
weights2 = ones(1,length(data2));
weights2(1,row_indices2) = 0;
row_indices12 = datasample(1:length(data1), cell_sample_size/2, 'Weights', weights1, 'Replace', false);
row_indices22 = datasample(1:length(data2), cell_sample_size/2, 'Weights', weights2, 'Replace', false);
row_indices_all2 = [row_indices12'; row_indices22'];
rand_data12 = data1(row_indices12,:);
rand_data22 = data2(row_indices22,:);
rand_data_all2 = [rand_data12; rand_data22];
%% Sample more different cells
weights1 = ones(1,length(data1));
weights1(1,[row_indices1 row_indices12]) = 0;
weights2 = ones(1,length(data2));
weights2(1,[row_indices2 row_indices22]) = 0;
row_indices1 = datasample(1:length(data1), cell_sample_size/2, 'Weights', weights1, 'Replace', false);
row_indices2 = datasample(1:length(data2), cell_sample_size/2, 'Weights', weights2, 'Replace', false);
row_indices = [row_indices1'; row_indices2'];
rand_data1 = data1(row_indices1,:);
rand_data2 = data2(row_indices2,:);
rand_data = [rand_data1; rand_data2];
%%
%% Sample initial conditions
format short g      % displays more decimals
init_sample_size = 100;           % change this if changing number of starting initial conditions
ref = [0.1775, log(2)/9, 0.5188, 1E-6, 0.004, 14, 0.08];     % reference point (from 2021 paper)
x0 = [0.1775, log(2)/9, 0.5188, 1E-6, 0.004];       % parameters with 0.1x-10x constraint
delay0 = 14;
nfkb_init0 = 0.08;
% Get random number [-1,1] as power of 10 to then multiply reference point
rand_x = 10.^(-1 + (1-(-1)).*rand(init_sample_size,length(x0))).*x0;
% rand_x = log10(rand_x);
% Custom constraints
rand_delays = 7 + (28-7).*rand(init_sample_size,1);
rand_nfkb_inits = 0.04 + (0.3-0.04).*rand(init_sample_size,1);
rand_inits = [rand_x rand_delays rand_nfkb_inits];
%%
%% Set up simulations
lb = [0.01775, 0.1*log(2)/9, 0.05188, 1E-7, 0.0004, 7, 0.04];     % lower bound
ub = [1.775, 10*log(2)/9, 5.188, 1E-5, 0.04, 28, 0.3];        % upper bound
options.Display = 'off';
options.PlotFcns = [];
options.MaxIter = 100;

params_fit = zeros(cell_sample_size, numel(ref), init_sample_size);
%%
%% Run the parameter fitting algorithm
first_half_c = cell_sample_size/2;
sec_half_c = cell_sample_size/2+1;
% mycluster = parcluster;
% mycluster.NumWorkers = 50;
% mycluster.parpool(50);
parfor i = 1:size(rand_inits,1)
 for c = 1:first_half_c
     xi = invert_constraint(rand_inits(i,:),lb,ub);
     fun = @(p) opt_rmsd(p, sheet_one, row_indices1(c)+1);
     [p_fit, fval, exitflag, output] = fminsearch(fun, xi, options);
     p_fit = constrain(p_fit,lb,ub);
     params_fit(c,:,i) = p_fit;
 end
 for c = sec_half_c:cell_sample_size
     xi = invert_constraint(rand_inits(i,:),lb,ub);
     fun = @(p) opt_rmsd(p, sheet_two, row_indices2(c-cell_sample_size/2)+1);
     [p_fit, fval, exitflag, output] = fminsearch(fun, xi, options);
     p_fit = constrain(p_fit,lb,ub);
     params_fit(c,:,i) = p_fit;
 end
 disp("parallel iteration " + string(i) + " completed.")
end
%%
%% Re-simulate with new parameters
end_rmsd = zeros(cell_sample_size, init_sample_size);
raw_rmsd = zeros(cell_sample_size, init_sample_size);
top_rmsds = zeros(cell_sample_size, 10);
best_params = zeros(cell_sample_size, numel(ref), 10);
mkdir(folder_name);
parfor c = 1:cell_sample_size
 sim_data = zeros(init_sample_size, 98);
 cur_end_rmsd = zeros(1,init_sample_size);
 for i = 1:init_sample_size
     p_mod = [6 3 params_fit(c,1,i); 15 1 params_fit(c,2,i); 16 1 params_fit(c,2,i); 66 1 params_fit(c,3,i); 68 1 params_fit(c,4,i); 75 1 params_fit(c,5,i); 6 4 params_fit(c,6,i)];
     init_mod = {'NFkB', params_fit(c,7,i)};
     doses = 100;
     dose_scale = 1/1500; % Convert to uM. Pam3CSK molecular weight: 1.5KDa
     names = {'NFkBn', 'IkBaNFkBn'};
     options = struct;
     options.DEBUG = 0;
     options.SIM_TIME = 8*60 + 5;
     % Simulate (only need to equilibrate on first iteration)
     output = [];
     if isempty(output)
         try
             [t,x,simdata2] = nfkbSimulate({'Pam3CSK',doses*dose_scale},names, p_mod, init_mod,options);
         catch
             warning('on');
             warning('nfkbSimulate failed. Proceeding to set rmsd to nan.');
             warning('off');
             cur_end_rmsd(i) = NaN;
             end_rmsd(c,i) = cur_end_rmsd(i);
             raw_rmsd(c,i) = cur_end_rmsd(i);
             continue
         end
     else
         options.STEADY_STATE = simdata2.STEADY_STATE;
         [~,x] = nfkbSimulate({'Pam3CSK',doses*dose_scale},names, p_mod, {},options);
     end
     output = cat(3,output,x);
     sim_nfkb_curve = squeeze(output(:,strcmp(names,'NFkBn'),:));
     sim_ikba_curve = squeeze(output(:,strcmp(names,'IkBaNFkBn'),:));
     sim_nfkb_curve = sim_nfkb_curve + sim_ikba_curve;
     sim_nfkb_curve = sim_nfkb_curve(1 : 5 : end);
     % bring sim to same zero point as experimental data
     zero = sim_nfkb_curve(1,1) - rand_data(c,1);
     sim_nfkb_curve = sim_nfkb_curve - zero;
     sim_data(i,:) = sim_nfkb_curve;
     dif = sim_nfkb_curve - rand_data(c,:)';
     raw_rmsd(c,i) = rms(dif, 'omitnan');
     dif(1:48,1) = dif(1:48,1).*sqrt(2);
     ratio = params_fit(c,:,i)./ref;
     ratio(1,1:5) = log10(ratio(1,1:5));
     ratio(1,6:7) = log2(ratio(1,6:7));
     cur_end_rmsd(i) = rms(dif, 'omitnan') + 0.1*rms(ratio);
     end_rmsd(c,i) = cur_end_rmsd(i);
 end
 % Get top 10% rmsds
 [cur_end_rmsd,indices] = sort(cur_end_rmsd);
 top_rmsds(c,:) = cur_end_rmsd(1:10);
 top_sim_data = sim_data(indices',:);
 top_sim_data = top_sim_data(1:10,:);
 top_params = params_fit(c,:,indices);
 top_params = top_params(1,:,1:10);
 best_params(c,:,:) = top_params;
 % Plot data
 fig = figure;
 sums = sum(sim_data,2,'omitnan');
 [~,I] = sort(sums, 'descend');
 ranked_sim = sim_data(I',:);
 plot(top_sim_data')
 hold on
 % Ax = gca;
 % Ax.XDisplayLabels = nan(size(Ax.XDisplayData));
 % Ax.YDisplayLabels = nan(size(Ax.YDisplayData));
 % colormap default
 % Plot experimental data
 plot(rand_data(c,:)', 'LineWIdth', 3)
 if c <= first_half_c
     plot_name = 'P3Ki0' + string(1) + '_row_#' + string(row_indices(c)+1);
 else
     plot_name = 'P3Ki0' + string(2) + '_row_#' + string(row_indices(c)+1);
 end
 print(fig, plot_name, '-dpng')
 plot_file = plot_name + '.png';
 movefile(plot_file, folder_name);
 hold off
end
% Plot data
% sums = sum(rand_data,2,'omitnan');
% [~,I] = sort(sums, 'descend');
% ranked_exp = rand_data(I',:);
% plot(rand_data')
% Ax = gca;
% Ax.XDisplayLabels = nan(size(Ax.XDisplayData));
% Ax.YDisplayLabels = nan(size(Ax.YDisplayData));
% Ax.Title = 'Experimental';
% colormap default
%%
%% Store relevant data
result = cell(1,7);      % (1) cells, (2) initial conditions, (3) param fits
result{1} = [row_indices1 + 1; row_indices2 + 1];
result{2} = rand_inits;
result{3} = params_fit;
result{4} = raw_rmsd;
result{5} = end_rmsd;
result{6} = top_rmsds;
result{7} = best_params;
save('result_cell', 'result');
tEnd = toc(tStart);
save('runtime', 'tEnd');
movefile('result_cell.mat', folder_name);
movefile('runtime.mat', folder_name);
