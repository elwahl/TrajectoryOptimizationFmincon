function [fval, bigX] = collocate(numStates, numControls, valuesAndBounds, t_0, t_f, timePoints, thrust, m_dot, normValues, costBalanceTimeImportance, attemptEllipticalOrbit, performRobustnessAnalysis_Thrust, performRobustnessAnalysis_Phi, performRobustnessAnalysis_FiringTime)
    %% Do some initialization of things we'll need.
    numNodes = size(timePoints, 2);
    numStuff = numStates + numControls;

    %% Create the bigX and lower and upper bound matrices that we'll need.
    bigX = zeros((numNodes * numStuff) + 1, 1);
    lBounds = zeros((numNodes * numStuff) + 1, 1);
    uBounds = zeros((numNodes * numStuff) + 1, 1);

    %% Iterate through the initial/final values and bounds that we were passed.
    for i = 0:(numStuff - 1)
        currentValuesAndBounds = valuesAndBounds(i + 1);
        % Do the initial values guess. -ELW
        if isfield(currentValuesAndBounds, 'overrideInitGuessValue') && isscalar(currentValuesAndBounds.overrideInitGuessValue)
            % This is for when we have an explicitly specified value we want for the initial guesses. -ELW
            bigX(((i * numNodes) + 1):((i + 1) * numNodes)) = (ones(1, numNodes) * currentValuesAndBounds.overrideInitGuessValue);
        else
            % This is for when we are assuming linear progression for the initial guesses. -ELW
            bigX(((i * numNodes) + 1):((i + 1) * numNodes)) = linspace((currentValuesAndBounds.initialLow + currentValuesAndBounds.initialHigh) / 2, (currentValuesAndBounds.finalLow + currentValuesAndBounds.finalHigh) / 2, numNodes);
        end
        % Populate the lower bounds array. -ELW
        lBounds((i * numNodes) + 1) = currentValuesAndBounds.initialLow;
        lBounds((i * numNodes) + 2 : ((i + 1) * numNodes) - 1) = currentValuesAndBounds.middleLow;
        lBounds((i + 1) * numNodes) = currentValuesAndBounds.finalLow;
        % Populate the upper bounds array. -ELW
        uBounds((i * numNodes) + 1) = currentValuesAndBounds.initialHigh;
        uBounds((i * numNodes) + 2 : ((i + 1) * numNodes) - 1) = currentValuesAndBounds.middleHigh;
        uBounds((i + 1) * numNodes) = currentValuesAndBounds.finalHigh;
    end

    % Special overrides! -ELW
    bigX(1 : floor(numNodes/2)) = linspace(0, 0.2, floor(numNodes/2));
    bigX(ceil(numNodes/2) : numNodes) = linspace(0.2, 0, numNodes - ceil(numNodes/2) + 1);
    if (attemptEllipticalOrbit)
        bigX(numNodes * 5) = 0;
        vDiffOffset = 1;
    else
        vDiffOffset = 0;
    end
    bigX((numNodes * (4 + vDiffOffset)) + 1 : (numNodes * (4 + vDiffOffset)) + floor(numNodes/2)) = linspace(0, pi, floor(numNodes/2));
    bigX((numNodes * (4 + vDiffOffset)) + ceil(numNodes/2) : (numNodes * (5 + vDiffOffset))) = linspace(-pi, 0, numNodes - ceil(numNodes/2) + 1);

    % Extra time-handling stuff. -ELW
    bigX(end) = (t_f + t_0) / 2;    % For free final time, t_0 and t_f represent the range of values we think t_f might actually take. -ELW
    lBounds(end) = 0;               % If t_f is free, then we should make the bounds be very loose. -ELW
    uBounds(end) = t_f * 2;         % We'll let this value be relatively large, but not completely out of context. -ELW

    %% See if we want to do any robustness analysis.
    % Set any values ahead of time so that we don't confuse fmincon too much by setting them in varying ways during program execution. -ELW
    randomization = zeros(numNodes, 3);
    if (performRobustnessAnalysis_Thrust)
        for i = 1 : numNodes
            %randomization(i, 1) = 1.001 - (rand * .002);    % 0.1 percent case. -ELW
            %randomization(i, 1) = 1.01 - (rand * .02);      % 1.0 percent case. -ELW
            %randomization(i, 1) = 1.05 - (rand * .1);       % 5.0 percent case. -ELW
            randomization(i, 1) = 1.1 - (rand * .2);        % 10.0 percent case. -ELW
        end
    end
    if (performRobustnessAnalysis_Phi)
        for i = 1 : numNodes
            %randomization(i, 2) = 1.001 - (rand * .002);    % 0.1 percent case. -ELW
            %randomization(i, 2) = 1.01 - (rand * .02);      % 1.0 percent case. -ELW
            %randomization(i, 2) = 1.05 - (rand * .1);       % 5.0 percent case. -ELW
            randomization(i, 2) = 1.1 - (rand * .2);        % 10.0 percent case. -ELW
        end
    end
    if (performRobustnessAnalysis_FiringTime)
        for i = 1 : numNodes
            %randomization(i, 3) = 1.001 - (rand * .002);    % 0.1 percent case. -ELW
            %randomization(i, 3) = 1.01 - (rand * .02);      % 1.0 percent case. -ELW
            %randomization(i, 3) = 1.05 - (rand * .1);       % 5.0 percent case. -ELW
            randomization(i, 3) = 1.1 - (rand * .2);        % 10.0 percent case. -ELW
        end
    end

    %% Set options, and actually run fmincon!
    tic;
    options = optimset('Algorithm','sqp','MaxFunEvals',100000000,'Display','iter-detailed','MaxIter',10000);
    [bigX, fval] = fmincon(@(X) objfun(X, numNodes, costBalanceTimeImportance), bigX, [], [], [], [], lBounds, uBounds, @(X) nonlconstr(X, numStates, numNodes, timePoints, thrust, m_dot, normValues, attemptEllipticalOrbit, performRobustnessAnalysis_Thrust, performRobustnessAnalysis_Phi, performRobustnessAnalysis_FiringTime, randomization), options);
    toc;        % Show how much time has passed. -ELW

    %% Display the results.
    missionNumberDays = bigX(end) * normValues.time / (24 * 3600);
    missionFuelConsumed = (bigX((3 * numNodes) + 1) - bigX(4 * numNodes)) * normValues.mass;
    disp(['Length of time for the mission: ', num2str(missionNumberDays), ' days']);
    disp(['Amount of fuel consumed for the mission: ', num2str(missionFuelConsumed), ' kg']);
    plotTimeSpan = linspace(0, missionNumberDays, 11);
    numSubPlots = 6 + vDiffOffset;
    currentSubPlot = 1;
    subplot(round(numSubPlots/2), 2, currentSubPlot), plot(bigX(1 : numNodes) * normValues.velocity,'r','LineWidth',2);
    grid on; xlabel('time (days)'); ylabel('km/s'); title('Radial Velocity');
    set(gca, 'XTickLabel', sprintf('%.2f|', plotTimeSpan));
    currentSubPlot = currentSubPlot + 1;
    subplot(round(numSubPlots/2), 2, currentSubPlot), plot(bigX(numNodes + 1 : 2 * numNodes) * normValues.velocity,'r','LineWidth',2);
    grid on; xlabel('time (days)'); ylabel('km/s'); title('Tangential Velocity');
    set(gca, 'XTickLabel', sprintf('%.2f|', plotTimeSpan));
    currentSubPlot = currentSubPlot + 1;
    subplot(round(numSubPlots/2), 2, currentSubPlot), plot(bigX(2 * numNodes + 1 : 3 * numNodes) * normValues.radius,'r','LineWidth',2);
    grid on; xlabel('time (days)'); ylabel('km'); title('Radius');
    set(gca, 'XTickLabel', sprintf('%.2f|', plotTimeSpan));
    currentSubPlot = currentSubPlot + 1;
    subplot(round(numSubPlots/2), 2, currentSubPlot), plot(bigX(3 * numNodes + 1 : 4 * numNodes) * normValues.mass,'r','LineWidth',2);
    grid on; xlabel('time (days)'); ylabel('kg'); title('Spacecraft Mass');
    set(gca, 'XTickLabel', sprintf('%.2f|', plotTimeSpan));
    currentSubPlot = currentSubPlot + 1;
    if (attemptEllipticalOrbit)
        subplot(round(numSubPlots/2), 2, currentSubPlot), plot(bigX(4 * numNodes + 1 : 5 * numNodes) * normValues.velocity,'r','LineWidth',2);
        grid on; xlabel('time (days)'); ylabel('km/s'); title('V diff');
        set(gca, 'XTickLabel', sprintf('%.2f|', plotTimeSpan));
        currentSubPlot = currentSubPlot + 1;
    end
    subplot(round(numSubPlots/2), 2, currentSubPlot), plot(bigX(numStates * numNodes + 1 : (numStates + 1) * numNodes),'r','LineWidth',2);
    grid on; xlabel('time (days)'); ylabel('radians'); title('phi');
    set(gca, 'XTickLabel', sprintf('%.2f|', plotTimeSpan));
    currentSubPlot = currentSubPlot + 1;
    subplot(round(numSubPlots/2), 2, currentSubPlot), plot(bigX((numStates + 1) * numNodes + 1 : (numStates + 2) * numNodes),'r','LineWidth',2);
    grid on; xlabel('time (days)'); ylabel('decimal percent'); title('Thruster Firing Percentage');
    set(gca, 'XTickLabel', sprintf('%.2f|', plotTimeSpan));

    if (performRobustnessAnalysis_Thrust || performRobustnessAnalysis_Phi || performRobustnessAnalysis_FiringTime)
        figure;
        subplot(3, 1, 1), plot(randomization(1 : numNodes, 1),'r','LineWidth',2);
        grid on; xlabel('time (days)'); ylabel('Decimal Percent'); title('Thrust Variation');
        set(gca, 'XTickLabel', sprintf('%.2f|', plotTimeSpan));
        subplot(3, 1, 2), plot(randomization(1 : numNodes, 2),'r','LineWidth',2);
        grid on; xlabel('time (days)'); ylabel('Decimal Percent'); title('Phi Variation');
        set(gca, 'XTickLabel', sprintf('%.2f|', plotTimeSpan));
        subplot(3, 1, 3), plot(randomization(1 : numNodes, 3),'r','LineWidth',2);
        grid on; xlabel('time (days)'); ylabel('Decimal Percent'); title('Firing Time Variation');
        set(gca, 'XTickLabel', sprintf('%.2f|', plotTimeSpan));
    end
end