function [F, Q, F_block, Q_block, EigVals, B_weights] = nipals_multiBlock(X_normalized,...
    num_PCs, tol, max_iter, deflationMethod)
% Function to implement the NIPALS method from Hanafi 2010
% Options to perform either Multiple Co-Inertia Analysis or Consensus Principle Component Analysis
% * Does NOT include any preprocessing steps on input data. 
% * Does NOT normalize output vectors.
% Inputs: 
%         * X_normalized    Vertical cell array of NORMALIZED data matrices (each has n rows)
%         * num_PCs         Number of PCs desired in output
%         * tol             Scalar tolerance in stopping criterion for NIPALS
%         * max_iter        Maximum number of iterations allowed in NIPALS algorithm
%         * deflationMethod Choice of block deflation method:
%                               'block' to use block loadings for deflation (MCIA - default)
%                               'global' to use global scores for deflation (CPCA)
% Outputs:
%         * F               Matrix of global scores
%         * Q               Matrix of global loadings
%         * F_block         Cell array of local score matrices  
%         * Q_block         Cell array of local loading matrices
%         * EigVal          Array of eigenvalues corresponding to first 'num_PCs' global scores
%         * B_weights       Matrix where i^th column has the block contributions to i^th global score


if nargin < 5
    deflationMethod = 'block';
end

% Getting variable and sample numbers for each dataset
num_datasets = length(X_normalized);
num_samples = size(X_normalized{1},1); % Extract number of samples in datasets (i.e. # rows) 
num_variables = zeros(1,num_datasets); % Array with number of variables in each block
for i = 1:num_datasets
    num_variables(1,i) = size(X_normalized{i},2);
end
total_vars = sum(num_variables,'all'); % Total number of variables across whole dataset

% Global score/loadings storage matrices
F = zeros(num_samples,num_PCs); % Matrix of global scores - each column is one global score

Q = zeros(total_vars,num_PCs); % Matrix of global loadings - each column is one global loading (F = X Q)

% Block score storage cell array 
F_block = cell(1,num_datasets); 
% CELL index specifies the dataset
% For cell index i: F_block{i} = [f_i^(1) | f_i^(2) | ... | f_i^(num_PCs)]
for i = 1:num_datasets
    F_block{i} = zeros(num_samples, num_PCs);
end

% Block loadings storage cell array 
Q_block = cell(1,num_datasets); 
% CELL index specifies the dataset
% For cell index i: Q_block{i} = [q_i^(1) | q_i^(2) | ... | q_i^(num_PCs)]  
for i = 1:num_datasets
    X_i = X_normalized{i};
    num_vars = size(X_i); num_vars = num_vars(2);
    Q_block{i} = zeros(num_vars, num_PCs);
end

%% NIPALS implementation (see Hanafi et. al. 2010)

X_deflated = X_normalized; % Creating copy of data to run deflation process

% Creating data super-matrix
X_super = X_normalized{1};
if num_datasets > 1
    for i = 2:num_datasets
        X_super = [X_super X_normalized{i}];
    end
end

% Creating eigenvalue list
EigVals = zeros(1,num_PCs);

% Creating block weight array
B_weights = zeros(num_datasets,num_PCs);

% Main NIPALS iteration
for j = 1:num_PCs % repeat deflation/iteration process for desired number of PCs

    %%% Initialization 
    T = zeros(num_samples,num_datasets); % matrix to store block scores
    
    block_loadings = cell(1,num_datasets); % cell array to store block loading vectors
    
    f = ones(num_samples,1); % arbitrary starting vector
   
    iter = 0; % iteration count
    as_old = 0; % stopping criterion variable
    err = 10; % tolerance metric
    if ~exist('max_iter')
        max_iter = 10000; % iteration limiter if not specified
    end
    
    %%% NIPALS iteration loop
    while err > tol
        for i = 1:num_datasets
            X_i = X_deflated{i}; % extract ith dataset
            
            q_i = X_i'*f/(f'*f); % compute block loadings 

            q_i = q_i/norm(q_i); % normalize loading vector

            f_i = X_i*q_i; % compute block scores

            T(:,i) = f_i; % add block score to matrix;
            
            block_loadings{i} = q_i; % store block loadings

        end

        w = T'*f/(f'*f); % compute global weights (relative to matrix of 

        w = w/norm(w); % normalize global weights

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
           fprintf('WARNING: exceeded maximum iteration threshold \n')
           break 
        end
    end 
    
    %%% Updating outputs not used in iteration
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
    
    % jth global score -->  jth column of matrix F
    F(:,j) = f;
    
    % jth set of block weights --> jth column of B_weights
    B_weights(:,j) = w; 
    
    % Computing global loadings
    Qi = Q_block{1};
    q_global = B_weights(1,j)*Qi(:,j);
    q_global= q_global';
    % compute global loadings by weighted concatenation of block loadings:
    if num_datasets > 1
        for i = 2:num_datasets
            Qi = Q_block{i};
            q_global = [q_global (B_weights(i,j)*Qi(:,j))';];
        end
    end
    % jth global loading --> jth column of matrix Q
    Q(:,j) = q_global;
    
    % Creating deflated data super-matrix
    X_super_deflated = X_deflated{1};
    if num_datasets > 1
        for i = 2:num_datasets
            X_super_deflated = [X_super_deflated X_deflated{i}];
        end
    end
     
    % Recording first eigenvalue of deflated super-matrix (corresponding covariance matrix)
    svs = svd(X_super_deflated);
    EigVals(j) = svs(1);
    
    fprintf(['Completed after ', num2str(iter), ' iterations \n']);
    
    %%% Deflation step
    % Deflate using either global scores (for CPCA) or block loadings (for MCIA)
    if strcmp('global',deflationMethod)
        % Deflation via global scores
        for i = 1:num_datasets
            X_i = X_deflated{i}; % extract ith dataset

            Xnew_i = X_i - (f*f'/(f'*f))*X_i;

            X_deflated{i} = Xnew_i;
        
        end
    else
        % Deflation via block loadings
        for i = 1:num_datasets
            X_i = X_deflated{i}; % extract ith dataset 
            q_i = Q_block{i}; % extract block loading matrix related to ith dataset

            Xnew_i = X_i - X_i*q_i(:,j)*(q_i(:,j)');

            X_deflated{i} = Xnew_i;
        end
    end
end


end