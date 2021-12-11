function [blocks_normalized, EigVals] = omicade_initialization(blocks_raw)
% Function to apply the preprocessing from the omicade4.R package for MCIA.
%   INPUTS:
%       * blocks_raw - cell array of data matrices for different omics types (in sample x variable format)
%                      *** Each matrix must contain same number of rows ***
%   OUTPUTS:
%       * blocks_normalized - cell array of data matrices with initialization applied


    %% Step 1: Making data non-negative
    blocks_normalized = blocks_raw;
    for i = 1:length(blocks_normalized)
        minN = min(min(blocks_normalized{i}));
        if(minN < 0)
           offset = floor(minN);
           blocks_normalized{i} = blocks_normalized{i} + abs(offset); 
        end
    end
    %% Performing row/column weightings 
    % Weights for later step:
    N = zeros(1,length(blocks_normalized)); % stores sums of all values in a block
    col_w = cell(length(blocks_normalized),1); % stores column sum contributions
    row_w = cell(length(blocks_normalized),1); % stores row sum contributions

    for i = 1:length(blocks_normalized)
        data = blocks_normalized{i}; 
        N(i) = sum(data,'all'); % NOTE: this is really poorly-conditioned for nearly-empty blocks
        if N(i) < 1e-3 
            warning(['Block ', num2str(i),' is nearly constant.'])
        end
        if N(i) == 0 
            error(['Block ', num2str(i),' is constant, cannot apply initialization.'])
        end
        col_w{i} = ones(1,size(data,1))*data/N(i); % column sums/total block sums
        row_w{i} = data*ones(size(data,2),1)/N(i); % row sums/total block sums

    end

    % Step 2: Divide each column by its sum
    for j = 1:length(blocks_normalized)
        data = blocks_normalized{j};
        for i = 1:size(data,2) 
            colsum = sum(data(:,i));
            if colsum < 1e-3 && colsum ~= 0 
                warning(['Block ', num2str(j),' column ',num2str(j),' is nearly zero.'])
            end
            if colsum ~= 0 
                data(:,i) = data(:,i).*(1/colsum);
            end
        end
        blocks_normalized{j} = data; 
    end

    % Step 3: Subtract sum contribution percent from each row 
    for i = 1:length(blocks_normalized)
        data = blocks_normalized{i};
        rowW = row_w{i};
        for j = 1:size(data,1)
            data(j,:) = data(j,:) - rowW(j);   
        end
        blocks_normalized{i} = data; 
    end


    % Step 4: multiply each block by number of rows 
    for i = 1:length(blocks_normalized)
        blocks_normalized{i} = blocks_normalized{i}*size(blocks_normalized{i},1); 
    end

    %% Computing block weights via singular values: 
    eval_list = cell(1); % array of singular values for normalized blocks (each block has a different row)

    for k = 1:length(blocks_normalized)
        df = blocks_normalized{k};

        % multiply each column by its weight:
        colweights = col_w{k};
        for i =1:size(df,2)
            df(:,i) = df(:,i).*sqrt(colweights(i));
        end
        
        % divide by square root of number of samples:
        df = df/sqrt(size(df,1));

        % Compute covariance matrix
        df = df*df';

        % Computing eigenvalues of covariance matrix (i.e. singular values) 
        [V, D] = eig(df);
        eigvals = diag(D);

        % Dropping eigenvalues/vectors less than tol
        tol = 1e-7;
        indices = find(abs(eigvals)<tol);
        for i =1:length(indices)
            j = indices(i);
            V(:,j) = [];
            eigvals(j) = []; % note: eigenvalues are ordered lowest to highest 
        end

        % Sorting eigenvalues highest to lowest 
        [eigvals, ind] = sort(eigvals,'descend');
        V = V(:, ind);
        EigVals = eigvals;
        
        % Saving eigvals:
        eval_list{k} = eigvals; 

    end

    %% Step 5: applying inertial block weighting

    % Divide each block by the square root of its total inertia
    % Inertia of block = sum of singular values
    for i = 1:length(blocks_normalized)
        blockW = 1/sum(eval_list{i});
        blocks_normalized{i} = blocks_normalized{i}*sqrt(blockW); 
    end

    %% Step 6: applying col and row weights

    % Weighting each table by sqrt( # of rows): 
    for k =1:length(blocks_normalized)
        blocks_normalized{k} = blocks_normalized{k}*1/sqrt(size(blocks_normalized{k},1));
    end 

    % Weighting columns:
    for k =1:length(blocks_normalized)
        df = blocks_normalized{k};
        colweights = col_w{k};
        for i =1:size(df,2)
            df(:,i) = df(:,i).*sqrt(colweights(i));
        end
        blocks_normalized{k} = df; 
    end 

end

