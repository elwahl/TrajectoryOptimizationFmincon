function f = objfun(bigX, numNodes, costBalanceTimeImportance)
    %% This is the cost function that we want to minimize.
    % Note that we can now do a balance between time and fuel consumption,
    % depending on what we set costBalanceTimeImportance to be, in the
    % range of 0 to 1. -ELW
    f = (costBalanceTimeImportance * bigX(end)) + ((1 - costBalanceTimeImportance) * (bigX((3 * numNodes) + 1) - bigX(4 * numNodes)));
end