function x = multiBlockPlot(global_coords,block_coords,clusIndices, clusNames)
% Plots a linked scatterplot from the first two global and block scores
% Intended to be used on joint embeddings from MCIA or CPCA, where underlying clusters are known
% NOTE: inputs may need normalization prior to plotting. 
% Inputs: 
%         * global_coords     nx2 matrix of first two global scores
%         * block_coords      cell array of nx2 matrices of first two block scores
%         * clusIndices       cell array of indices for each cluster
%                               {clus1_indices, clus2_indices, ...}
%         * clusNames         cell array of strings of cluster names
%                               {clus1_name, clus2_name,...} 

    x = [];% empty output

    blockMarkerList = {'o','s','^','*','x','+'};
    colorList = {'g','b','r','m','y','c'};
    numClusters = length(clusNames);
    numBlocks = max(size(block_coords));
    
    if numClusters > length(colorList)
       error('Number of clusters exceeds number of plotting colors') 
    end
    
    if numBlocks > length(blockMarkerList)
       error('Number of blocks exceeds number of plotting symbols') 
    end
    
    % Plotting global coordinates:
    figure()
    for i = 1:numClusters
        clusX_global = global_coords(clusIndices{i},1);
        clusY_global = global_coords(clusIndices{i},2);
        symb = [colorList{i},'.'];
        scatter(clusX_global,clusY_global,symb); hold on
        
        % Plotting block coordinates connected by a line
        for j = 1:numBlocks
            block_score = block_coords{j};
            clusX_block = block_score(clusIndices{i},1);
            clusY_block = block_score(clusIndices{i},2);
            for k = 1:length(clusX_block)
                plot([clusX_block(k) clusX_global(k)], [clusY_block(k) clusY_global(k)],...
                    ['-',colorList{i}],'HandleVisibility','off')
                plot(clusX_block(k),clusY_block(k),[colorList{i},blockMarkerList{j}],...
                    'HandleVisibility','off')
            end
        end
    end
    legend(clusNames,'Location','southeast')
    
    title('Plot of Projections onto First Two Scores')
    grid on;
end

