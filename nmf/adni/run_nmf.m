
% Script to run Aris' OPNMF implementation on the
% tau ROI data from ADNI.

% SET PATHS TO THINGS HERE
ARIS_CODE_PATH = '../brainparts-master';
INPUT_CSV_PATH = '../../derivatives/adni/nmf_matrix.csv';
OUTPUTDIR = 'results';

% SET NMF PARAMETERS HERE

Ks = 2:1:20;

w0 = '';
initMeth = '';
max_iter = '';
tol = '';
iter0 = '';
save_step = 1000000; % doesn't save intermediate results!
outputdir = OUTPUTDIR;

% import code
addpath(ARIS_CODE_PATH);

% load data
X = csvread(INPUT_CSV_PATH);

% run!

mkdir(OUTPUTDIR);

MATDIR = fullfile(OUTPUTDIR, 'mat');
mkdir(MATDIR);

recErrors = zeros(size(Ks));
for i = 1:length(Ks)

    disp(i);

    % run NMF
    K = Ks(i);
    [W, H] = opnmf(X, K, w0, initMeth, max_iter, tol, iter0, save_step, outputdir);
    
    % add normalized component projection
    Wnorm = W ./ sum(W, 1);
    Hnorm = X' * Wnorm;

    savepath = fullfile(MATDIR, sprintf('NumBases%d.mat', K));
    save(savepath, 'W', 'H', 'Wnorm', 'Hnorm');

    % Save recon error
    Xhat = W * H;
    E = (X - Xhat);
    recErrors(i) = norm(E, 'fro');

end
