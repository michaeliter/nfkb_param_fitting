function D = cell_pair_dist_calc(x, y, min_vals, max_vals, bin_widths, means, stds)
%x = rand(10, 7);
%y = rand(10, 7);

JSD = zeros(1, 7);
for i = 1:7
    %fx = ksdensity(x(:, i),pts(:, i));
    %fy = ksdensity(y(:, i),pts(:, i));
    %M = (fx+fy)/2;
    
    x(:, i) = (x(:, i) - means(i))/stds(i);
    y(:, i) = (y(:, i) - means(i))/stds(i);
    bins = min_vals(i):bin_widths(i):max_vals(i);
    
    fx = histcounts(x(:, i), bins, 'Normalization', 'pdf')*bin_widths(i);
    fy = histcounts(y(:, i), bins, 'Normalization', 'pdf')*bin_widths(i);
    M = (fx+fy)/2;
    
    select = fx > 0;
    fx = fx(select);
    M1 = M(select);
    
    KLD1 = sum(fx.*log(fx./M1));
    
    select = fy > 0;
    fy = fy(select);
    M2 = M(select);
    
    KLD2 = sum(fy.*log(fy./M2));
    
    JSD(i) = sqrt(0.5*(KLD1 + KLD2));
end

D = mean(JSD);
end
