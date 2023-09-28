% JSD 
% Preprocessing
load params_3d.mat
all_params = params_3d;
%calculate bin width with Freedman Diaconis rule
min_vals = zeros(1, 7);
max_vals = zeros(1, 7);
bin_widths = zeros(1, 7);
means = zeros(1, 7);
stds = zeros(1, 7);
numCells = size(cell_matrix, 1);

for i = 1:7
    all = all_params(:, i, :);
    all = all(:);
    means(i) = mean(all);
    stds(i) = std(all);
    all = (all - means(i))/stds(i);
    min_vals(i) = min(all);
    max_vals(i) = max(all);
    bin_width = 2*iqr(all)*((length(all)/numCells)^(-1/3));
    bins = linspace(min_vals(i), max_vals(i), round((max_vals(i)-min_vals(i))/bin_width)+1);
    bin_widths(i) = bins(2) - bins(1);
    %how many bins on average, might be too few to detect meaningful
    %differences (divide bin_widths by n)
end
% Calculate distances per cell
jsd_matrix = zeros(size(cell_matrix,1));

for c1 = 1:numCells
    for c2 = c1:numCells
        X = squeeze(cell_matrix(c1,:,:))';
        Y = squeeze(cell_matrix(c2,:,:))';
        try
            jsd_matrix(c1,c2) = cell_pair_dist(X, Y, min_vals, max_vals, bin_widths, means, stds);
        catch
            error("cell " + string(c1) + ", cell " + string(c2))
        end
    end
    if mod(c1,100) == 0
        disp("cell " + string(c1) + " done")
    end
end


jsd_matrix;
jsd_matrix2 = transpose(jsd_matrix);
jsd_matrix = jsd_matrix + jsd_matrix2