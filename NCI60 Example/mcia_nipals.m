% Function to implement MCIA via the NIPALS method 
% Inputs: 
%         * X_list      Vertical cell array of data matrices (n rows)
%         * num_PCs     Number of PCs desired in output
%         * tol         Scalar tolerance in stopping criterion for PCs
%         * normalize   Boolean if algorithm should apply normalization 
% Outputs:
%         * F           Matrix of global scores
%         * Q           Matrix of global loadings
%         * F_block     Cell array of local score matrices  
%         * Q_block     Cell array of local loading matrices



function [F, Q, F_block, Q_block] = mcia_nipals(X_list, num_PCs, tol, normalize)
num_datasets = length(X_list);

n = size(X_list{1});
n = n(1);  % Extract number of samples in datasets (i.e. # rows) 

% Global score/loadings storage matrices
F = zeros(n,num_PCs); % Matrix of global scores - each column is one global score

% Q = zeros(total_vars,num_PCs); % Matrix of global loadings - each column is one global loading (F = X Q)

% Block score storage cell array 
F_block = cell(1,num_datasets); 
% CELL index specifies the dataset
% For cell index i: F_block{i} = [f_i^(1) | f_i^(2) | ... | f_i^(num_PCs)]
for i = 1:num_datasets
    F_block{i} = zeros(n, num_PCs);
end

% Block loadings storage cell array 
Q_block = cell(1,num_datasets); 
% CELL index specifies the dataset
% For cell index i: Q_block{i} = [q_i^(1) | q_i^(2) | ... | q_i^(num_PCs)]  
for i = 1:num_datasets
    X_i = X_list{i};
    num_vars = size(X_i); num_vars = num_vars(2);
    Q_block{i} = zeros(num_vars, num_PCs);
end


%% Normalization

if normalize
    X_normalized = cell(1); % Cell array containing normalized datasets.
    
    % Variable-level normalization:
    for i =1:length(X_list)
        X = X_list{i};
        for j = 1:length(X)
            X(:,j) = X(:,j)- mean(X(:,j)); % unit column mean        
            X(:,j) = X(:,j)/sqrt(var(X(:,j))); % unit column variance
        end
        X_normalized{i} = X;
    end
    
    % Block-level normalization
    % NOTE: division by sqrt(length) proposed in Tenenhaus & Tenenhaus, 2014 
    for i =1:length(X_normalized)
        X_normalized{i} = X_normalized{i}/sqrt(length(X_normalized{i}));
    end
else
    X_normalized = X_list; % No normalization applied
end


%% NIPALS implementation (see Hanafi et. al. 2010)

X_deflated = X_normalized; % Creating copy of data to run deflation process

for j = 1:num_PCs % repeat deflation/iteration process for desired number of PCs

    %%% Initialization 
    T = zeros(n,num_datasets); % matrix to store block scores
    
    block_loadings = cell(1,num_datasets); % cell array to store block loading vectors
    
    f = ones(n,1); % arbitrary starting vector
   
    iter = 0; % iteration count
    as_old = 0; % stopping criterion variable
    err = 10; % tolerance metric
    max_iter = 1000; % iteration limiter 

    %%% Iteration
    while err > tol
        for i = 1:num_datasets
            X_i = X_deflated{i}; % extract ith dataset
            
            q_i = X_i'*f/(f'*f); % compute block loadings 

            q_i = q_i/norm(q_i); % normalize loading vector

            f_i = X_i*q_i; % compute block scores

            T(:,i) = f_i; % add block score to matrix;
            
            block_loadings{i} = q_i; % store block loadings

        end

        w = T'*f/(f'*f); % compute global WEIGHTS (not loadings)

        w = w/norm(w); % normalize global WEIGHTS

        f = T*w; % update global score 

        % calculating stopping criterion:
        fnorm = f/sqrt(cov(f)); % unit variance global scores

        % iterating stopping metric
        as = 0; 
        for i = 1:num_datasets
            C = cov(T(:,i),fnorm);
            C = C(2)^2; % taking the covariance between qnorm and block scores. 
            as = as + C;
        end

        % checking stopping metric vs. stopping criterion
        err = abs(as - as_old); % error metric
        as_old = as;

        iter = iter+1;
        
        if iter > max_iter
           fprintf('Warning, exceeded maximum iteration threshold \n')
           break 
        end
    end 
    
    % Updating outputs:
    for i =1:num_datasets
        % record jth order block score in matrix relating to block i:
        Fi = F_block{i}; 
        Fi(:,j) = T(:,i);
        F_block{i} = Fi;
        
        % record jth order block loading in matrix relating to block i:
        Qi = Q_block{i};
        Qi(:,j) = block_loadings{i}; 
        Q_block{i} = Qi;
    end
    
    F(:,j) = fnorm; % ith global score => ith column of matrix F 
    
    fprintf(['Completed after ', num2str(iter), ' iterations \n']);
    
    %%% Deflation step
    % Deflate using block loadings as in Hassani et. al. 2013
    
    for i = 1:num_datasets
        X_i = X_deflated{i}; % extract ith dataset
        q_i = Q_block{i}; % extract block loading matrix related to ith dataset
        
        Xnew_i = X_i - X_i*q_i(:,j)*(q_i(:,j)');
        
        X_deflated{i} = Xnew_i;
        
    end
        
end
    

% DEBUG - activate global loadings matrices later. 
Q = 0;

end