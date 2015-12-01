%% Clear out any old stuff that's still laying around. -ELW
clear all;
close all;
clc;

% This defines the cost balance of time vs. fuel consumption.
% 1 = we only care about time, 0 = we only care about fuel, anywhere inbetween will be a blance of the two. -ELW
%costBalanceTimeImportance = 1;
costBalanceTimeImportance = 0.5;
%costBalanceTimeImportance = 0;

% Do we want to attempt optimization of an elliptical orbit (vs. the default circular orbit)? 1 = yes, 0 = no -ELW
% Solving for an elliptical orbit does not really work at this time, so better just to leave it at 0! -ELW
attemptEllipticalOrbit = 0;

% Do we want to do some robustness analysis by introducing perturbations? 1 = yes, 0 = no -ELW
performRobustnessAnalysis_Thrust = 0;
performRobustnessAnalysis_Phi = 0;
performRobustnessAnalysis_FiringTime = 0;

numTimeSlices = 199;                         % The number of time slices. -ELW
%numTimeSlices = 99;                         % The number of time slices. -ELW
%numTimeSlices = 4;                         % The number of time slices. -ELW

%% Put some initial/final and boundary values in for the different states/controls that we're going to use.
% This includes:
% V_r_i (the radial velocity of the spacecraft)
% V_theta_i (the tangential velocity of the spacecraft)
% r_i (the radius of the spacecraft)
% m_i (the mass of the spacecraft)
% V_diff_i (the difference between Mars' velocity at r_i and the spacecraft's velocity at r_i)
% phi_i (the direction of the spacecraft's thrust)
% delta_t_i (the percentage of each time slice that the thruster is firing)
% tau (the total transit time from Earth to Mars)

if (attemptEllipticalOrbit)
    numStates = 5;
else
    numStates = 4;
end
numControls = 2;        % Does not include tau. -ELW

%% Set some other values we'll need.
earthVelocityRadial = 0;                    % Radial velocity of the Earth, assuming circular orbit. (km/s) -ELW
earthVelocityTangential = 29.784736;        % Tangential velocity of the Earth, assuming circular orbit. (km/s) -ELW
earthRadius = 149597870;                    % Radius of the Earth's orbit, assuming circular orbit. (km) -ELW
marsRadiusMin = 206700000;                  % Minimum radius of Mars' orbit. (km) -ELW
marsVelocityMax = 26.5;                     % Maximum possible velocity of Mars. (km/s) -ELW
if (attemptEllipticalOrbit)
    marsRadiusMax = 249200000;              % Maximum radius of Mars' orbit (elliptical). (km) -ELW
    marsVelocityMin = 21.95;                % Minimum possible velocity of Mars (elliptical). (km/s) -ELW
else
    marsRadiusMax = marsRadiusMin;          % Maximum radius of Mars' orbit (circular). (km) -ELW
    marsVelocityMin = marsVelocityMax;      % Minimum possible velocity of Mars (circular). (km/s) -ELW
end
marsTransitionTimeMin = 16070400;           % Minimum amount of time we think the transition to Mars might take. (s) -ELW
marsTransitionTimeMax = 24105600;           % Maximum amount of time we think the transition to Mars might take. (s) -ELW
muSun = 132712440018;                       % Gravitational parameter for the sun. (km^3/s^2) -ELW
epsilonMars = -291.1638;                    % Mechanical energy of Mars. (km^2/s^2) -ELW

% These values will change depending on the type of engine used. -ELW
spacecraftMassDry = 3295;                   % Mass of the spacecraft without any fuel. (Bryson/Ho example, kg) -ELW
spacecraftMassWet = 4000;                   % Mass of the spacecraft when fully fueled. (Bryson/Ho example engine - kg) -ELW
spacecraftThrust = 0.003781;                % Thrust of the spacecraft. (Bryson/Ho example, kN) -ELW
spacecraftM_dot = 0.000067866;              % M_dot of the spacecraft. (Bryson/Ho example, kg/s) -ELW

%spacecraftMassDry = 2000;                   % Mass of the spacecraft without any fuel. (Electric engine, kg) -ELW
%spacecraftMassWet = 2600;                   % Mass of the spacecraft when fully fueled. (Electric engine, continuous thrust - kg) -ELW
%spacecraftThrust = 0.0005;                  % Thrust of the spacecraft. (Electric engine, kN) -ELW
%spacecraftM_dot = 0.000017;                 % M_dot of the spacecraft. (Electric engine, kg/s) -ELW

%spacecraftMassDry = 5000;                   % Mass of the spacecraft without any fuel. (Chemical engine, kg) -ELW
%spacecraftMassWet = 1364000;                % Mass of the spacecraft when fully fueled. (Chemical engine, max 10 min thrust - kg) -ELW
%spacecraftThrust = 10000;                   % Thrust of the spacecraft. (Chemical engine, kN) -ELW
%spacecraftM_dot = 2265;                     % M_dot of the spacecraft. (Chemical engine, kg/s) -ELW

%% Set the values we'll use for different normalizations.
normValues.time = (marsTransitionTimeMin + marsTransitionTimeMax) / 2;
normValues.radius = earthRadius;
normValues.gravity = .00980665;             % Gravitational value at Earth's surface. (km/s^2) -ELW
normValues.velocity = earthVelocityTangential;
normValues.mass = spacecraftMassWet;
normValues.thrust = normValues.mass * normValues.gravity;
normValues.mass_dot = normValues.gravity / normValues.velocity;
normValues.muSun = muSun / ((normValues.radius ^ 2) * normValues.gravity);
normValues.epsilonMars = epsilonMars / (normValues.velocity ^ 2);

%% Now we take the previous values we defined and normalize them all.
earthVelocityRadial = earthVelocityRadial / normValues.velocity;
earthVelocityTangential = earthVelocityTangential / normValues.velocity;
earthRadius = earthRadius / normValues.radius;
marsRadiusMin = marsRadiusMin / normValues.radius;
marsRadiusMax = marsRadiusMax / normValues.radius;
marsVelocityMin = marsVelocityMin / normValues.velocity;
marsVelocityMax = marsVelocityMax / normValues.velocity;
marsTransitionTimeMin = marsTransitionTimeMin / normValues.time;
marsTransitionTimeMax = marsTransitionTimeMax / normValues.time;
spacecraftMassDry = spacecraftMassDry / normValues.mass;
spacecraftMassWet = spacecraftMassWet / normValues.mass;
spacecraftThrust = spacecraftThrust / normValues.thrust;
spacecraftM_dot = spacecraftM_dot / normValues.mass_dot;

% Conditions for V_r_i:
valuesAndBounds.initialLow = earthVelocityRadial;
valuesAndBounds.initialHigh = earthVelocityRadial;
valuesAndBounds.middleLow = -10;
valuesAndBounds.middleHigh = 10;
%valuesAndBounds.middleLow = -3;
%valuesAndBounds.middleHigh = 3;
if (attemptEllipticalOrbit)
    valuesAndBounds.finalLow = -0.25 * marsVelocityMax;
    valuesAndBounds.finalHigh = 0.25 * marsVelocityMax;
else
    %valuesAndBounds.finalLow = -0.001;
    %valuesAndBounds.finalHigh = 0.001;
    valuesAndBounds.finalLow = 0;
    valuesAndBounds.finalHigh = 0;
end

% Conditions for V_theta_i:
valuesAndBounds(2).initialLow = earthVelocityTangential;
valuesAndBounds(2).initialHigh = earthVelocityTangential;
valuesAndBounds(2).middleLow = -10;
valuesAndBounds(2).middleHigh = 10;
%valuesAndBounds(2).middleLow = -3;
%valuesAndBounds(2).middleHigh = 3;
%valuesAndBounds(2).finalLow = marsVelocityMin - 0.001;
%valuesAndBounds(2).finalHigh = marsVelocityMax + 0.001;
valuesAndBounds(2).finalLow = marsVelocityMin;
valuesAndBounds(2).finalHigh = marsVelocityMax;

% Conditions for r_i:
valuesAndBounds(3).initialLow = earthRadius;
valuesAndBounds(3).initialHigh = earthRadius;
valuesAndBounds(3).middleLow = earthRadius;
%valuesAndBounds(3).middleHigh = 1.01 * marsRadiusMax;
valuesAndBounds(3).middleHigh = marsRadiusMax;
valuesAndBounds(3).finalLow = marsRadiusMin;
valuesAndBounds(3).finalHigh = marsRadiusMax;

% Conditions for m_i:
valuesAndBounds(4).initialLow = spacecraftMassWet;
valuesAndBounds(4).initialHigh = spacecraftMassWet;
valuesAndBounds(4).middleLow = spacecraftMassDry;
valuesAndBounds(4).middleHigh = spacecraftMassWet;
valuesAndBounds(4).finalLow = spacecraftMassDry;
valuesAndBounds(4).finalHigh = spacecraftMassWet;

if (attemptEllipticalOrbit)
    % Conditions for V_diff_i:
    valuesAndBounds(5).initialLow = -10;
    valuesAndBounds(5).initialHigh = 10;
    valuesAndBounds(5).middleLow = -10;
    valuesAndBounds(5).middleHigh = 10;
    valuesAndBounds(5).finalLow = -0.001;
    valuesAndBounds(5).finalHigh = 0.001;
    valuesAndBounds(5).overrideInitGuessValue = 1;
end

% Conditions for phi_i:
valuesAndBounds(numStates + 1).initialLow = -pi;
valuesAndBounds(numStates + 1).initialHigh = pi;
valuesAndBounds(numStates + 1).middleLow = -pi;
valuesAndBounds(numStates + 1).middleHigh = pi;
valuesAndBounds(numStates + 1).finalLow = -pi;
valuesAndBounds(numStates + 1).finalHigh = pi;
valuesAndBounds(numStates + 1).overrideInitGuessValue = pi / 2;

% Conditions for delta_t_i:
valuesAndBounds(numStates + 2).initialLow = 0;
valuesAndBounds(numStates + 2).initialHigh = 1;
valuesAndBounds(numStates + 2).middleLow = 0;
valuesAndBounds(numStates + 2).middleHigh = 1;
valuesAndBounds(numStates + 2).finalLow = 0;
valuesAndBounds(numStates + 2).finalHigh = 1;
%valuesAndBounds(numStates + 2).overrideInitGuessValue = 0.5;
valuesAndBounds(numStates + 2).overrideInitGuessValue = 0.75;

% Since we're handling the free end time case, all time slices should run from 0 to 1,
% as they will be scaled later. Note that the progression of values does NOT have to be linear. -ELW
timePoints = linspace(0, 1, numTimeSlices + 1);

%% Call the main collocate function. -ELW
minVal = collocate(numStates, numControls, valuesAndBounds, marsTransitionTimeMin, marsTransitionTimeMax, timePoints, spacecraftThrust, spacecraftM_dot, normValues, costBalanceTimeImportance, attemptEllipticalOrbit, performRobustnessAnalysis_Thrust, performRobustnessAnalysis_Phi, performRobustnessAnalysis_FiringTime);
