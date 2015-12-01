function [C, Ceq] = nonlconstr_robust(bigX, numStates, numNodes, timePoints, thrust, m_dot, normValues, attemptEllipticalOrbit, performRobustnessAnalysis_Thrust, performRobustnessAnalysis_Phi, performRobustnessAnalysis_FiringTime, randomization)
    % This function is a special version of nonlconstr for robustness
    % analysis. 
    % This function introduces upto 5% uncertainty into control phi and
    % thruster firing time. 
    % By introducing uncertainty into control parameters we can see if
    % mission can be carried out when things are not perfect!
   
    %% Nonlinear inequality constraints
    C = [];

    %% Nonlinear equality constraints - the functions we evaluate should want to be equal to 0.
    Ceq = zeros(numNodes, numStates);
    tau = bigX(end);            % Factor we want to scale the time slice by. (Effectively, the total time.) -ELW

    % Handle the special case of the first iteration. -ELW
    V_r_i = bigX(1);                                % (Guessed) value of V_r_i. -ELW
    V_theta_i = bigX(numNodes + 1);                 % (Guessed) value of V_theta_i. -ELW
    r_i = bigX((numNodes * 2) + 1);                 % (Guessed) value of r_i. -ELW
    m_i = bigX((numNodes * 3) + 1);                 % (Guessed) value of m_i. -ELW

    % The initial values are whatever fmincon sets them to, subject to the
    % lbound and ubound values. Therefore, Ceq for initial values should
    % always just be 0. -ELW
    %Ceq(1, 1) = V_r_i;                              % Take care of initial V_r_i. -ELW
    %Ceq(1, 2) = V_theta_i - 1;                      % Take care of initial V_theta_i. -ELW
    %Ceq(1, 3) = r_i - 1;                            % Take care of r_i. -ELW
    %Ceq(1, 4) = m_i - 1;                            % Take care of m_i. -ELW

    if (attemptEllipticalOrbit)
        % Take care of initial V_diff_i. -ELW
        V_diff_i = bigX((numNodes * 4) + 1);            % (Guessed) value of V_diff_i. -ELW
        Ceq(1, 5) = V_diff_i - sqrt(abs(2 * ((1 / r_i) + normValues.epsilonMars))) + sqrt((V_r_i ^ 2) + (V_theta_i ^ 2));
    end

    % Now we do the more general iterations. -ELW
    thrustOriginal = thrust;
    loopOffset = 1;         % 0 or 1 depeding on which way we want to calculate things. -ELW
    for n = (loopOffset + 1) : (numNodes - 1 + loopOffset)
        indexN = n + 1 - loopOffset;
        indexNm1 = n - loopOffset;

        % First get some basic values we'll need to do the rest of the calculations. -ELW
        deltaTau = abs((timePoints(indexN) - timePoints(indexNm1)) * tau);       % timePoints array contains percentage values for slices of tau. -ELW
        t = abs(bigX((numNodes * (numStates + 1)) + indexNm1) * deltaTau);     % t is the amount of time for the thruster to fire. -ELW
        V_r_im1 = bigX(indexNm1);                           % Value of V_r_i-1. -ELW
        V_theta_im1 = bigX(numNodes + indexNm1);            % Value of V_theta_i-1. -ELW
        r_im1 = bigX((numNodes * 2) + indexNm1);            % Value of r_i-1. -ELW
        m_im1 = bigX((numNodes * 3) + indexNm1);            % Value of m_i-1. -ELW
        phi_im1 = bigX((numNodes * numStates) + indexNm1);  % Value of phi_i-1. -ELW

        % Vary any of our calculation values by some random, pre-determined value, if desired. -ELW
        if (performRobustnessAnalysis_Thrust)
            thrust = thrustOriginal * randomization(n, 1);
        end
        if (performRobustnessAnalysis_Phi)
            phi_im1 = phi_im1 * randomization(n, 2);
        end
        if (performRobustnessAnalysis_FiringTime)
            t = t * randomization(n, 3);
            t = min(1, t);      % t cannot go outside of the 0 to 1 range. -ELW
        end

        % Now calculate the rate change values, first for the time period where the thruster is firing. -ELW
        m_new = m_im1 - (m_dot * t);
        V_r_i_dot = ((V_theta_im1 ^ 2) * normValues.muSun / r_im1) - (normValues.muSun / r_im1^2) + ((thrust * sin(phi_im1))/ m_new);
        V_theta_i_dot = ((-1 * V_r_im1 * V_theta_im1 * normValues.muSun) / r_im1) + ((thrust * cos(phi_im1))/ m_new);

        % Find the intermediate values. -ELW
        V_r_i_calc = V_r_im1 + (V_r_i_dot * t * normValues.time * normValues.gravity / normValues.velocity);
        V_theta_i_calc = V_theta_im1 + (V_theta_i_dot * t * normValues.time * normValues.gravity / normValues.velocity);
        r_i_calc = r_im1 + (V_r_i_calc * t * normValues.time * normValues.velocity / normValues.radius);

        % Then calculate the change values for when the thruster is not firing. -ELW
        V_r_i_dot2 = ((V_theta_i_calc ^ 2) * normValues.muSun / r_i_calc) - (normValues.muSun / r_i_calc^2);
        V_theta_i_dot2 = ((-1 * V_r_i_calc * V_theta_i_calc * normValues.muSun) / r_i_calc);

        % And then some final calculated values. -ELW
        V_r_i_calc2 = V_r_i_calc + (V_r_i_dot2 * (deltaTau - t) * normValues.time * normValues.gravity / normValues.velocity);
        V_theta_i_calc2 = V_theta_i_calc + (V_theta_i_dot2 * (deltaTau - t) * normValues.time * normValues.gravity / normValues.velocity);
        r_i_calc2 = r_i_calc + (V_r_i_calc2 * (deltaTau - t) * normValues.time * normValues.velocity / normValues.radius);

        % Last, do the state change equalities. -ELW
        V_r_i = bigX(indexN);                                % (Guessed) value of V_r_i. -ELW
        V_theta_i = bigX(numNodes + indexN);                 % (Guessed) value of V_theta_i. -ELW
        r_i = bigX((numNodes * 2) + indexN);                 % (Guessed) value of r_i. -ELW
        m_i = bigX((numNodes * 3) + indexN);                 % (Guessed) value of m_i. -ELW

        % The first equation takes care of V_r_i. -ELW
        Ceq(indexN, 1) = V_r_i - V_r_i_calc2;
        % The second takes care of V_theta_i. -ELW
        Ceq(indexN, 2) = V_theta_i - V_theta_i_calc2;
        % The third takes care of r_i. -ELW
        Ceq(indexN, 3) = r_i - r_i_calc2;
        % The fourth takes care of m_i. -ELW
        Ceq(indexN, 4) = m_i - m_new;
        if (attemptEllipticalOrbit)
            % The fifth takes care of V_diff_i. -ELW
            V_diff_i = bigX((numNodes * 4) + indexN);            % (Guessed) value of V_diff_i. -ELW
            Ceq(indexN, 5) = V_diff_i - sqrt(abs(2 *((1 / r_i) + normValues.epsilonMars))) + sqrt((V_r_i_calc2 ^ 2) + (V_theta_i_calc2 ^ 2));
        end
    end
end

