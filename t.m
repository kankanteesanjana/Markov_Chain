clc;
clear all;
close all;

% Load data (uphill and downhill sessions)
load('uphillSession.mat');
load('downhillSession.mat');

% Extract speed and slope data
uphillSpeeds = uphillSessions(:, 4);
uphillSlopes = uphillSessions(:, 7);
downhillSpeeds = downhillSessions(:, 4);
downhillSlopes = downhillSessions(:, 7);

% Define number of bins for speed
numSpeedBins = 20;

% Set the range for valid data
validRange = [0 7.6]; 

% Filter out values outside the valid range for uphill speeds
filteredUphillSpeeds = uphillSpeeds(uphillSpeeds >= validRange(1) & uphillSpeeds <= validRange(2));

% Filter out values outside the valid range for downhill speeds
filteredDownhillSpeeds = downhillSpeeds(downhillSpeeds >= validRange(1) & downhillSpeeds <= validRange(2));

% Classify slopes and count occurrences for uphill data using custom tolerance
tolerance = 1e-5;

% Initialize for unique slopes and their counts for uphill data
uniqueUpSlopes = uphillSlopes(1);  % Initialize with the first slope
slopeCountsUp = 1;  % Initialize count for the first slope

for i = 2:length(uphillSlopes)
    % Check if the current slope is not within the tolerance of any previous slopes
    if all(abs(uphillSlopes(i) - uniqueUpSlopes) > tolerance)
        uniqueUpSlopes = [uniqueUpSlopes; uphillSlopes(i)];
        slopeCountsUp = [slopeCountsUp; 1];  % Initialize count for the new slope
    else
        % If the slope is already in uniqueUpSlopes, increment its count
        index = find(abs(uniqueUpSlopes - uphillSlopes(i)) < tolerance);
        slopeCountsUp(index) = slopeCountsUp(index) + 1;
    end
end

% Classify slopes and count occurrences for downhill data using custom tolerance
tolerance = 1e-5;

% Initialize containers for unique slopes and their counts for downhill data
uniqueDownhillSlopes = downhillSlopes(1);  % Initialize with the first slope
slopeCountsDown = 1;  % Initialize count for the first slope

for i = 2:length(downhillSlopes)
    % Check if the current slope is not within the tolerance of any previous slopes
    if all(abs(downhillSlopes(i) - uniqueDownhillSlopes) > tolerance)
        uniqueDownhillSlopes = [uniqueDownhillSlopes; downhillSlopes(i)];
        slopeCountsDown = [slopeCountsDown; 1];  % Initialize count for the new slope
    else
        % If the slope is already in uniqueDownhillSlopes, increment its count
        index = find(abs(uniqueDownhillSlopes - downhillSlopes(i)) < tolerance);
        slopeCountsDown(index) = slopeCountsDown(index) + 1;
    end
end

% Store results in variables
uphillResults = [uniqueUpSlopes, slopeCountsUp];
downhillResults = [uniqueDownhillSlopes, slopeCountsDown];

% Discretize speeds into states
speedBins = linspace(min([filteredUphillSpeeds; filteredDownhillSpeeds]), ...
    max([filteredUphillSpeeds; filteredDownhillSpeeds]), numSpeedBins + 1);

uphillSpeedStates = discretize(filteredUphillSpeeds, speedBins);
downhillSpeedStates = discretize(filteredDownhillSpeeds, speedBins);

% Initialize a structure to store slope and corresponding speeds for uphill
slopeSpeedInfoUp = struct('slope', {}, 'originalIndices', {});

% Iterate over unique slope values for uphill
for i = 1:length(uniqueUpSlopes)
    % Extract the current unique slope for uphill
    slope = uniqueUpSlopes(i);

    % Find indices of occurrences of this slope in the original uphill data
    indices = find(abs(uphillSlopes - slope) < tolerance);

    % Check if indices are not empty before storing information
    if ~isempty(indices)
        % Store information in the structure for uphill
        slopeSpeedInfoUp(end + 1).slope = slope;
        slopeSpeedInfoUp(end).originalIndices = indices;
    end
end


% Initialize a structure to store slope and corresponding speeds for downhill
slopeSpeedInfoDown = struct('slope', {}, 'originalIndices', {});

% Iterate over unique slope values for downhill
for i = 1:length(uniqueDownhillSlopes)
    % Extract the current unique slope for downhill
    slope = uniqueDownhillSlopes(i);

    % Find indices of occurrences of this slope in the original downhill data
    indices = find(abs(downhillSlopes - slope) < tolerance);

    % Check if indices are not empty before storing information
    if ~isempty(indices)
        % Store information in the structure for downhill
        slopeSpeedInfoDown(end + 1).slope = slope;
        slopeSpeedInfoDown(end).originalIndices = indices;
    end
end

% Generate mappings for uphill and downhill
uphillMapping = createMapping(uphillSpeedStates, slopeSpeedInfoUp);
downhillMapping = createMapping(downhillSpeedStates, slopeSpeedInfoDown);

numUniqueSlopesUp = length(uniqueUpSlopes);
numUniqueSlopesDown = length(uniqueDownhillSlopes);

% Generate transition matrices for uphill
uphillTransitionMatrices = cell(1, numUniqueSlopesUp);
for i = 1:numUniqueSlopesUp
    originalIndices = uphillMapping(i).originalIndices;
    states = uphillMapping(i).states;

    % Pass the correct number of states to the function
    uphillTransitionMatrices{i} = generateTransitionMatrix(states, numSpeedBins);
end

% Generate transition matrices for downhill
downhillTransitionMatrices = cell(1, numUniqueSlopesDown);
for i = 1:numUniqueSlopesDown
    originalIndices = downhillMapping(i).originalIndices;
    states = downhillMapping(i).states;

    % Pass the correct number of states to the function
    downhillTransitionMatrices{i} = generateTransitionMatrix(states, numSpeedBins);
end

totalNumberOfSteps = 500;  

% Simulate Markov chain for uphill
[simulatedMarkovChainUphill, simulatedSpeedsBySlopeUphill] = simulateMarkovChainUphill(uphillTransitionMatrices, slopeSpeedInfoUp, filteredUphillSpeeds, totalNumberOfSteps, numSpeedBins);


% Simulate Markov chain for downhill
[simulatedMarkovChainDownhill, simulatedSpeedsBySlopeDownhill] = simulateMarkovChainDownhill(downhillTransitionMatrices, slopeSpeedInfoDown, filteredDownhillSpeeds, totalNumberOfSteps, numSpeedBins);


% Calculate the stationary distribution using eigenvector method
T_Uphill = uphillTransitionMatrices{2}; 
[V_Uphill, D_Uphill] = eigs(T_Uphill', 1, 'largestabs');
%[V_Uphill, D_Uphill] = eigs(T_Uphill', size(T_Uphill, 1));

% Extract the eigenvector corresponding to the eigenvalue 1
stationaryDistribution_Uphill = V_Uphill;

% Normalize the stationary distribution
stationaryDistribution_Uphill = stationaryDistribution_Uphill / sum(stationaryDistribution_Uphill);

% Calculate probabilities after a long simulation
numBins = 100; 
edges = linspace(min(simulatedMarkovChainUphill), max(simulatedMarkovChainUphill), numBins + 1);
binIndices = discretize(simulatedMarkovChainUphill, edges);

% Use histcounts to calculate probabilities
finalStateFreq_Uphill = histcounts(simulatedMarkovChainUphill, edges) / numel(simulatedMarkovChainUphill);

% Normalize the final state frequencies
finalStateFreq_Uphill = finalStateFreq_Uphill / sum(finalStateFreq_Uphill);

% Check if the probabilities are approximately the same for Uphill
tolerance_Uphill = 1e-5;
areProbabilitiesApproxSame_Uphill = all(abs(stationaryDistribution_Uphill - finalStateFreq_Uphill) < tolerance_Uphill);

if areProbabilitiesApproxSame_Uphill
    disp('The stationary distribution and probabilities are approximately the same for Uphill.');
else
    disp('The stationary distribution and probabilities are not the same for Uphill.');
end

 % For Uphill
 figure;
 subplot(2,1,1);
 plot(uphillSpeeds(1:min(totalNumberOfSteps, length(uphillSpeeds))), 'LineWidth', 2);
 title('Original Speeds - Uphill');
 xlabel('Time Steps');
 ylabel('Speed');
 
 subplot(2,1,2);
 plot(simulatedMarkovChainUphill, 'LineWidth', 2);
 title('Simulated Markov Chain - Uphill');
 xlabel('Time Steps');
 ylabel('Speed');
 
 % For Downhill
 figure;
 subplot(2,1,1);
 plot(downhillSpeeds(1:min(totalNumberOfSteps, length(downhillSpeeds))), 'LineWidth', 2);
 title('Original Speeds - Downhill');
 xlabel('Time Steps');
 ylabel('Speed');

 subplot(2,1,2);
 plot(simulatedMarkovChainDownhill, 'LineWidth', 2);
 title('Simulated Markov Chain - Downhill');
 xlabel('Time Steps');
 ylabel('Speed');
 
 
% % Plot results for uphill simulation
% figure('Name', 'Simulated Speeds and Transition Matrices for Uphill');
% 
% % Plot simulated speeds
% subplot(2, 1, 1);
% plot(simulatedMarkovChainUphill, 'LineWidth', 2);
% title('Simulated Speeds for Uphill');
% xlabel('Time Steps');
% ylabel('Speed');
% 
% % Display transition matrices for uphill
% subplot(2, 1, 2);
% uphillConcatenatedMatrix = cat(1, uphillTransitionMatrices{:});
% imagesc(uphillConcatenatedMatrix);
% colorbar;
% title('Transition Matrices for Uphill');
% xlabel('Current State');
% ylabel('Next State');
% 
% % Plot results for downhill simulation
% figure('Name', 'Simulated Speeds and Transition Matrices for Downhill');
% 
% % Plot simulated speeds
% subplot(2, 1, 1);
% plot(simulatedMarkovChainDownhill, 'LineWidth', 2);
% title('Simulated Speeds for Downhill');
% xlabel('Time Steps');
% ylabel('Speed');
% 
% % Display transition matrices for downhill
% subplot(2, 1, 2);
% downhillConcatenatedMatrix = cat(1, downhillTransitionMatrices{:});
% imagesc(downhillConcatenatedMatrix);
% colorbar;
% title('Transition Matrices for Downhill');
% xlabel('Current State');
% ylabel('Next State');

% Filtered Uphill Speeds
filteredUphillSpeeds = uphillSpeeds(uphillSpeeds >= validRange(1) & uphillSpeeds <= validRange(2));

% Filtered Downhill Speeds
filteredDownhillSpeeds = downhillSpeeds(downhillSpeeds >= validRange(1) & downhillSpeeds <= validRange(2));

% For Uphill
figure;
subplot(3,1,1);
plot(filteredUphillSpeeds(1:min(totalNumberOfSteps, length(filteredUphillSpeeds))), 'LineWidth', 2);
title('Filtered Uphill Speeds');
xlabel('Time Steps');
ylabel('Speed');

subplot(3,1,2);
plot(simulatedMarkovChainUphill, 'LineWidth', 2);
title('Simulated Markov Chain - Uphill');
xlabel('Time Steps');
ylabel('Speed');

subplot(3,1,3);
plot(uphillSpeeds(1:min(totalNumberOfSteps, length(uphillSpeeds))), 'LineWidth', 2);
title('Original Speeds - Uphill');
xlabel('Time Steps');
ylabel('Speed');

% For Downhill
figure;
subplot(3,1,1);
plot(filteredDownhillSpeeds(1:min(totalNumberOfSteps, length(filteredDownhillSpeeds))), 'LineWidth', 2);
title('Filtered Downhill Speeds');
xlabel('Time Steps');
ylabel('Speed');

subplot(3,1,2);
plot(simulatedMarkovChainDownhill, 'LineWidth', 2);
title('Simulated Markov Chain - Downhill');
xlabel('Time Steps');
ylabel('Speed');

subplot(3,1,3);
plot(downhillSpeeds(1:min(totalNumberOfSteps, length(downhillSpeeds))), 'LineWidth', 2);
title('Original Speeds - downhill');
xlabel('Time Steps');
ylabel('Speed');




function mapping = createMapping(speedStates, slopeSpeedInfo)
    mapping = struct('slope', {}, 'originalIndices', {}, 'states', {});

    for i = 1:length(slopeSpeedInfo)
        slope = slopeSpeedInfo(i).slope;
        originalIndices = slopeSpeedInfo(i).originalIndices;

        % Check if originalIndices are within the bounds of speedStates
        validIndices = originalIndices(originalIndices <= length(speedStates));

        % Extract the states corresponding to the valid indices
        states = speedStates(validIndices);

        % Store information in the structure
        mapping(end + 1).slope = slope;
        mapping(end).originalIndices = validIndices;
        mapping(end).states = states;
    end
end

function transitionMatrix = generateTransitionMatrix(states, numStates)
    % Initialize transition matrix with zeros
    transitionMatrix = zeros(numStates);

    % Count occurrences of transitions
    for i = 1:length(states) - 1
        current = states(i);
        next = states(i + 1);
        transitionMatrix(current, next) = transitionMatrix(current, next) + 1;
    end

    % Normalize transition matrix
    rowSum = sum(transitionMatrix, 2);

    % Find non-zero elements (transitions that occurred)
    nonZeroElements = transitionMatrix > 0;

    % Normalize only the non-zero elements
    for row = 1:numStates
        if rowSum(row) ~= 0
            transitionMatrix(row, nonZeroElements(row, :)) = transitionMatrix(row, nonZeroElements(row, :)) / rowSum(row);
        end
    end

    % Handle NaN values 
    transitionMatrix(isnan(transitionMatrix)) = 0;

    % Handle zero sums 
    zeroSums = rowSum == 0;
    transitionMatrix(zeroSums, :) = 0;
end

