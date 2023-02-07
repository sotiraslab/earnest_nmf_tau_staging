
% Script to run Aris' OPNMF implementation on the
% tau ROI data from ADNI.

% SET PATHS TO THINGS HERE
ARIS_CODE_PATH = '../brainparts-master';
FULL_X_PATH = '../../derivatives/oasis3/nmf_matrix.csv';
SPLITS_PATH = '../../derivatives/oasis3/reproducibility_split_indices.csv';
OUTPUTDIR = 'results';

% SET NMF PARAMETERS HERE

Ks = 2:1:20;

w0 = '';
initMeth = '';
max_iter = '';
tol = '';
iter0 = '';
save_step = 1000000; % doesn't save intermediate results!
outputdir = '';

% import code
addpath(ARIS_CODE_PATH);

% load data
X = csvread(FULL_X_PATH);
SPLITS = csvread(SPLITS_PATH);
SHAPE = size(SPLITS);
N = SHAPE(2);

%%

% run!
REP_DIR = fullfile(OUTPUTDIR, 'repeat_reproducibility');

for n = 1:N
    ITERATION_DIR = fullfile(REP_DIR, sprintf('repeat%d', n));
    ITERATION_A = fullfile(ITERATION_DIR, 'repA');
    ITERATION_B = fullfile(ITERATION_DIR, 'repB'); 
    mkdir(ITERATION_DIR);
    mkdir(ITERATION_A);
    mkdir(ITERATION_B);

    idx = SPLITS(:, n);
    X_A = X(:, idx == 1);
    X_B = X(:, idx == 2);

    for i = 1:length(Ks)

        disp(i);
    
        % run NMF
        K = Ks(i);
    
        % rep A
        [W, H] = opnmf(X_A, K, w0, initMeth, max_iter, tol, iter0, save_step, outputdir);
        savepath = fullfile(ITERATION_A, sprintf('NumBases%d.mat', K));
        save(savepath, 'W', 'H');
    
        % rep B
        [W, H] = opnmf(X_B, K, w0, initMeth, max_iter, tol, iter0, save_step, outputdir);
        savepath = fullfile(ITERATION_B, sprintf('NumBases%d.mat', K));
        save(savepath, 'W', 'H');
    end
            
end




