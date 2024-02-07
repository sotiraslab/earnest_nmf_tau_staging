
% Modified version of Aris' code to run reproducibility analyses
% tweaked for my setup, running on ROIs

% See original code here:
%  https://github.com/sotiraslab/aris_nmf_analyses/blob/main/evaluateReproducibility.m

% Parameters - dependencies
PATH_ARIS_CODE = '../../scripts';

addpath(PATH_ARIS_CODE);

% Parameters - I/O
pathDir1 = '../../nmf/adni/results/mat';
pathDir2 = '../../nmf/oasis3/results/mat';
outputdir = 'adni_v_oasis_matching';

% Main function

% get list of results files from the first split
% note that this will not be naturally sorted
listing = dir(fullfile(pathDir1, '*.mat'));
filenames = {listing(:).name};

% extract the rank from each results, and get their natural order
basenums = regexp(filenames, '\d+', 'match'); % this is a nested cell array
basenums = horzcat(basenums{:}); % this is a non-nested cell array
basenums = cellfun(@str2num, basenums); % this is a numerical vector
[sortedBasisNum, idx] = sort(basenums);
sortedFileNames = filenames(idx);

% get full file names
fullfiles1 = fullfile(pathDir1, sortedFileNames);
fullfiles2 = fullfile(pathDir2, sortedFileNames);

numDifBases = numel(filenames);
ARI=zeros(numDifBases,1);

% create output
mkdir(outputdir);
matchdir = fullfile(outputdir, 'matching');
mkdir(matchdir);

for exp=1:numDifBases      
    p1 = fullfiles1{exp};
    p2 = fullfiles2{exp};
    disp(sprintf('Comparing files: \n  %s\n  %s', p1, p2));

    resSplit1 = load(p1) ;
    resSplit2 = load(p2) ;
    
    W1 = resSplit1.Wnorm;
    W2 = resSplit2.Wnorm;
    
    % calculate inner products
    inner_product = W1'*W2 ;
    
    % take a distance
    dist = 2*(1 - inner_product) ;
    
    % find correspondences
    [Matching,~] = Hungarian(dist);
    [~,idx_hug1]=max(Matching,[],2);

    % save correspondences
    outfile = fullfile(matchdir, sprintf('Match%d.mat', sortedBasisNum(exp)));
    save(outfile, 'Matching', 'idx_hug1');
    
    % overlap - hungarian
    overlap{exp} = zeros(sortedBasisNum(exp),1) ;
    for b=1:sortedBasisNum(exp)
        overlap{exp}(b) = inner_product(b,idx_hug1(b));
    end
    
    % overlap with best
    overlap_best{exp} = max(inner_product,[],2) ;
    
    % also evaluate overlap based on adjusted Rand Index    
    [~,clustering1] = max(W1,[],2);
    [~,clustering2] = max(W2,[],2);
    ARI(exp) = adjusted_rand_index(clustering1,clustering2);
    
end

meanInner=cellfun(@(x) mean(x),overlap,'UniformOutput', false);
medianInner=cellfun(@(x) median(x),overlap,'UniformOutput', false);

% save results
mkdir(outputdir)
save(fullfile(outputdir, 'reproducibilityResults.mat'), 'meanInner','medianInner','ARI','sortedBasisNum');
