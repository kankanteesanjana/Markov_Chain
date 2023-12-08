function [simulatedMarkovChain, simulatedSpeedsBySlope, stationaryDistribution] = simulateMarkovChainUphill(transitionMatrices, slopeSpeedInfo, filteredSpeeds, totalNumberOfSteps, numSpeedBins)
    numSlopes = length(slopeSpeedInfo);

    % Initialize a cell array to store simulated speeds for each slope
    simulatedSpeedsBySlope = cell(1, numSlopes);

    % Initialize an array to store simulated speeds
    simulatedMarkovChain = zeros(1, totalNumberOfSteps);

    % Function to choose the next speed index based on the transition probabilities
    chooseNextSpeedIndex = @(transitionProbabilities, numSpeeds, currentSpeedIndex, validIndices) chooseNextSpeedIndexFunction(transitionProbabilities, numSpeeds, currentSpeedIndex, validIndices);

    % Generate simulated speeds for each slope
    for i = 1:numSlopes
        slope = slopeSpeedInfo(i).slope;
        occurrences = length(slopeSpeedInfo(i).originalIndices);

        % Initialize a cell array to store simulated speeds for each slope
        simulatedSpeedsBySlope{i} = simulateChain(transitionMatrices{i}, occurrences, filteredSpeeds, numSpeedBins);

        % Update simulated Markov chain based on the simulated indices
        simulatedIndices = slopeSpeedInfo(i).originalIndices;

        % Ensure the simulated indices are within the total number of steps
        validIndices = simulatedIndices <= totalNumberOfSteps & simulatedIndices >= 1;

        % Update the simulated Markov chain only for valid indices
        simulatedMarkovChain(simulatedIndices(validIndices)) = simulatedSpeedsBySlope{i}(validIndices);
    end

end


function nextIndex = chooseNextSpeedIndexFunction(transitionProbabilities, numSpeeds, currentSpeedIndex, validIndices)
    cumulativeProbabilities = cumsum(transitionProbabilities);
    randomValue = rand();

    % Ensure the generated random index is within the bounds
    if all(cumulativeProbabilities == 0)
        nextIndex = randi([1, numSpeeds]);
    else
        % Choose the next index using the cumulative probabilities
        nextIndex = find(cumulativeProbabilities >= randomValue, 1);

        % Apply constraint: Ensure the next index is within the valid indices
        nextIndex = intersect(nextIndex, validIndices);

        % If no valid index is found, choose randomly from valid indices
        if isempty(nextIndex)
            nextIndex = randsample(validIndices, 1);
        end
    end
end


function simulatedSpeeds = simulateChain(transitionMatrix, occurrences, filteredSpeeds, numSpeedBins)
    numSpeeds = size(transitionMatrix, 1);

    % Initialize the simulated speeds array
    simulatedSpeeds = zeros(1, occurrences);

    % Choose a random starting speed index
    currentSpeedIndex = randi([1, numSpeeds]);

    for step = 1:occurrences
        % Store the current speed in the simulated speeds array
        simulatedSpeeds(step) = filteredSpeeds(currentSpeedIndex);

        % Choose the next speed index based on the transition probabilities
        transitionProbabilities = transitionMatrix(currentSpeedIndex, :);

        % Ensure transition probabilities are not empty
        if ~isempty(transitionProbabilities)
            % Get valid indices based on some criteria (e.g., nearby states)
            validIndices = getValidIndices(currentSpeedIndex, numSpeedBins);

            % Choose the next index with the constraint on valid indices
            nextIndex = chooseNextSpeedIndexFunction(transitionProbabilities, numSpeeds, currentSpeedIndex, validIndices);

            % Update the current speed index
            currentSpeedIndex = nextIndex;
        else
            % Handle the case of empty transition probabilities
            break;
        end
    end
end

function validIndices = getValidIndices(currentSpeedIndex, numSpeedBins)
    % Consider nearby states within a range of 5 indices to the left and right
    range = 100;
    
    % Ensure that the range is within bounds
    leftBound = max(1, currentSpeedIndex - range);
    rightBound = min(numSpeedBins, currentSpeedIndex + range);

    % Generate valid indices
    validIndices = leftBound:rightBound;
end
