%% generate a matrix of fake neural data
% what we want is a matrix that is Conditions x Neurons x Time
% these are simulated "condition-averaged" responses

% define our dimensionality 
T = 500;
C = 4;
N = 50;
% matrix to hold responses
response_matrix = zeros( C, N, T );

% indexing vectors
t = 1:T;
conditions = 1:C;

% group the neurons

% randomize the ordering of neurons
rand_numbers = rand(N, 1);
[~, neuron_order] = sort(rand_numbers);

%

groups = { neuron_order([1:25]), ...
           neuron_order([26:50]) };


% define firing pattern and condition dependence for first group
base_pattern{1} = sin( 2*pi*t / T);
condition_multiples{1} = [-2 1 0 4];

% define firing pattern and condition dependence for secondgroup
base_pattern{2} = cos( 2*pi*t / T);
condition_multiples{2} = [1 5 -7 2];


%% make the data
% use the info above to make fake neural responses
for group_number = 1:2
    % define the base response per condition for this group
    for nc = conditions
        responses{ group_number }( nc, : ) = base_pattern{ group_number } ...
            * condition_multiples{ group_number }( nc );
    end

    % each individual neuron's response is a random multiple of the base response
    for ineuron = 1:numel( groups{ group_number } )
        neuron = groups{ group_number }( ineuron );
        neuron_multiple = 1;%randn(1);

        for nc = conditions
            response_matrix( nc, neuron, : ) = responses{ group_number }( nc, : )...
                * neuron_multiple;
        end
    end
end

%% use the Tensor class to analyze this data
t = Tensor.Tensor;
t.t = response_matrix;
t.t = t.removeMarginalMeans();

subplot(1,2,1)
cmatrices = t.calcCovarianceMatrices();
imagesc( cmatrices.cy );
colorbar
title('unsorted matrix');


% get an ordering for the neuron correlation matrix
order = Utils.getCrossCorrOrdering( cmatrices.cy, 0.5 );
subplot(1,2,2)
imagesc( cmatrices.cy( order, order ) )
colorbar
title('sorted matrix');
