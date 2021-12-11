% Script to implement MCIA on a cut version of the NCI-60 cancer cell line dataset
% Example data taken from Omicade4 (Meng et. al. 2014)
% Dataset has three different omics types:
%       * mrna.csv  - 
%       * miRNA.csv - 
%       * prot.csv  - 

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
% The first two block scores are plotted in different shapes, connected to global scores
% With omicade initialization and 'block' deflation, replicates output of 'omicade' package

global_xcoords = GS_norm(:,1); 
global_ycoords = GS_norm(:,2);

CNSx = global_xcoords(1:6);
CNSy = global_ycoords(1:6);

LEUx = global_xcoords(7:12);
LEUy = global_ycoords(7:12);

MEx = global_xcoords(13:21);
MEy = global_ycoords(13:21);

datablockMarkers = {'o','s','^'}; 

clf; figure(1);
scatter(CNSx,CNSy, 'g.');
% Flipping axes if necessary:
set(gca, 'YDir','reverse'); % Flip y axis
% set(gca, 'XDir','reverse'); % Flip x axis
hold on;
scatter(LEUx,LEUy,'b.');
scatter(MEx,MEy,'r.');
legend('CNS','Leukaemia','Melanoma','Location','southeast')

for i =1:num_blocks
    block_score = BS_norm{i};
    block_x = block_score(:,1);
    block_y = block_score(:,2);
    
    CNS_x = block_x(1:6);
    CNS_y = block_y(1:6);
    for j = 1:length(CNS_x)
       plot([CNS_x(j) CNSx(j)], [CNS_y(j) CNSy(j)],'-g','HandleVisibility','off')
       plot(CNS_x(j),CNS_y(j),['g',datablockMarkers{i}],'HandleVisibility','off')
    end
    
    LEU_x = block_x(7:12);
    LEU_y = block_y(7:12);
    for j = 1:length(LEU_x)
       plot([LEU_x(j) LEUx(j)], [LEU_y(j) LEUy(j)],'-b','HandleVisibility','off')
       plot(LEU_x(j),LEU_y(j),['b',datablockMarkers{i}],'HandleVisibility','off')
    end
    
    ME_x = block_x(13:21);
    ME_y = block_y(13:21);
    for j = 1:length(ME_x)
       plot([ME_x(j) MEx(j)], [ME_y(j) MEy(j)],'-r','HandleVisibility','off')
       plot(ME_x(j),ME_y(j),['r',datablockMarkers{i}],'HandleVisibility','off')
    end
end
title('Plot of Projections onto First Two Scores')
grid on;

%% Plotting singular value decline
figure()
bar(evals.^2)
title('Singular Value Decline')
