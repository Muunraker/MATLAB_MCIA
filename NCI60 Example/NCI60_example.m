% Script to implement MCIA on a cut version of the NCI-60 cancer cell line dataset
% Example data taken from Omicade4 (Meng et. al. 2014)
% Dataset has three different omics types:
%       * mrna
%       * miRNA
%       * proteins

addpath('..\Functions\') % Path to MCIA functions
dataPath = '.\Data\'; % Path to folder containing data

%% Importing data:
redo_import = 1;
if redo_import
    mrna = readtable([dataPath,'mrna.csv'],'ReadRowNames',true);
    miRNA = readtable([dataPath,'miRNA.csv'],'ReadRowNames',true);
    prot = readtable([dataPath,'prot.csv'],'ReadRowNames',true);
end
    
% Creating cell array containing data:
blocks_raw = cell(1);
blocks_raw{1} = mrna{:,:}'; % Transpose ensures each dataset has same number of rows.
blocks_raw{2} = miRNA{:,:}';
blocks_raw{3} = prot{:,:}';


%% Data Preprocessing
% Pick either the preprocessing from omicade4 (Meng et. al 2014) or RGCCA (Tennenhaus & Tennenhaus, 2014)
% Omicade4 preprocessing:
blocks_normalized = omicade_initialization(blocks_raw);

% RGCCA preprocessing:
% blocks_normalized = RGCCA_initialization(blocks_raw);

%% Performing dimensionality-reduction via NIPALS
num_PCs = 10; % number of embeddings desired 
[Global_scores, Global_loadings, Block_scores, Block_loadings,evals,B_weights] = ...
    nipals_multiBlock(blocks_normalized,num_PCs,1e-14,10000, 'block');

%% Normalizing output
num_blocks = length(blocks_raw);
% Enforcing unit variance in global scores
GS_norm = Global_scores;
GL_norm = Global_loadings; 
for i = 1:num_PCs
    GS_norm(:,i) = GS_norm(:,i)/sqrt(var(GS_norm(:,i)));
    GL_norm(:,i) = GL_norm(:,i)/sqrt(var(GL_norm(:,i)));
end

% Enforcing unit variance in block scores
BS_norm = Block_scores;
BL_norm = Block_loadings;
for i = 1:num_blocks
    block_score = BS_norm{i};
    block_loading = BL_norm{i};
    for j = 1:num_PCs
        block_score(:,j) = block_score(:,j)/sqrt(var(block_score(:,j)));
        block_loading(:,j) = block_loading(:,j)/sqrt(var(block_loading(:,j)));
    end
    BS_norm{i} = block_score;
    BL_norm{i} = block_loading;
end 


%% Plotting Results
% Plotting first two global scores (principal components) in sample space
% Each block score is plotted with a different shape, connected by a line to the global score

% Compare to Figure 2A in (Meng et. al., 2016) if using 'omicade_initialization' 
% and 'block' deflation.

% Extracting first two global and block scores
global_coords = GS_norm(:,1:2);
block_coords = cell(1,num_blocks);
for i =1:num_blocks
    block_coords{i} = BS_norm{i}(:,1:2);
end

% Known clusters for cell lines in the NCI60 dataset: CNS, Leukemia, and Melanoma
CNS = 1:6; LEU = 7:12; ME = 13:21; 
clusIndices = {CNS,LEU,ME};
clusNames = {'CNS','Leukaemia','Melanoma'};

multiBlockPlot(global_coords,block_coords,clusIndices,clusNames);
% set(gca, 'YDir','reverse'); Fliping y axis to match plot in (Meng et. al., 2016) 

%% Plotting singular value decline at each iteration
figure()
bar(evals.^2)
title('Singular Value Decline')
