function ari = adjusted_rand_index(x,y)

    % Function for calculating the Adjusted Rand Index
    % based on the formulation here:
    % https://en.wikipedia.org/wiki/Rand_index#Adjusted_Rand_index

    % this is basically a translation of Python code,
    % so may not be implemented optimally for MATLAB

    % x - array of integer cluster labels
    % y - array of integer cluster labels
    
    if length(x) ~= length(y)
        error('Inputs x and y must have equal length');
    end
    
    n = length(x);
    clusters_x = unique(x);
    clusters_y = unique(y);
    cmat = zeros(length(clusters_x), length(clusters_y));
    
    for i=1:length(clusters_x)
        for j=1:length(clusters_y)
            xlabel = clusters_x(i);
            ylabel = clusters_y(j);
            cmat(i, j) = sum((x == xlabel) & (y == ylabel));
        end
    end
    
    a = sum(cmat, 2);
    b = sum(cmat, 1);
    
    cmat_comb = sum(arrayfun(@(val) mynchoosek(val, 2), cmat), 'all');
    a_comb = sum(arrayfun(@(val) mynchoosek(val, 2), a));
    b_comb = sum(arrayfun(@(val) mynchoosek(val, 2), b));
    n_comb = nchoosek(n, 2);
    
    top = cmat_comb - ((a_comb * b_comb)/n_comb);
    bot = (0.5 * (a_comb + b_comb)) - ((a_comb * b_comb)/n_comb);
    
    ari = top/bot;
    
    function ans = mynchoosek(n, k)
        if (k == 0) | (n < k)
            ans = 0;
        else
            ans = nchoosek(n, k);
        end
    end
end

