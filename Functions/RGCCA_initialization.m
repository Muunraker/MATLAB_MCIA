function [blocks_normalized] = RGCCA_initialization(blocks_raw)
% Function to implement preprocessing suggested for RGCCA (Tenenhaus & Tenenhaus, 2014).
%   INPUTS:
%       * blocks_raw - cell array of data matrices for different omics types (in sample x variable format)
%                      *** Each matrix must contain same number of rows ***
%   OUTPUTS:
%       * blocks_normalized - cell array of data matrices with initialization applied

blocks_normalized = cell(1); % Cell array containing normalized datasets.

% Variable-level normalization:
for i =1:length(blocks_raw)
    X = blocks_raw{i};
    for j = 1:length(X)
        X(:,j) = X(:,j)- mean(X(:,j)); % subtract unit column mean        
        X(:,j) = X(:,j)/sqrt(var(X(:,j))); % divide by sqrt(unit column variance)
    end
    blocks_normalized{i} = X;
end

% Block-level normalization
% NOTE: division by sqrt(length) proposed in Tenenhaus & Tenenhaus, 2014 
for i =1:length(blocks_normalized)
    blocks_normalized{i} = blocks_normalized{i}/sqrt(length(blocks_normalized{i}));
end

end

