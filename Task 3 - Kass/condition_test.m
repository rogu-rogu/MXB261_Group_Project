function sample_success = condition_test(samples, m, X0, T)
% CONDITION_TEST - Tests an array of 3-tuple parameter values and outputs
% the proportion of successful trials (for conditions provided by the
% assignment).
%
%   sample_success = condition_test(samples,m,X0,T) iterates through the
%   rows of samples, running m many trials for each parameter 3-tuple and
%   determining success if cond1 or cond2 is true.
%   * cond1: S(T) ∈ [S(0) - 20, S(10) - 10]
%   * cond2: S(T) < S(0) - 100
%
%   Input Arguments
%     samples - (n x 3) matrix, elements ∈[0,1] 
%       3-tuple parameter samples, with all elements ranging from 0 to 1.
%     m - integer
%       Number of trials; the number of times the Gillespie stochastic
%       simulation algorithm generates a new set of values.
%     X0 - (1 x 4) vector
%       The initial condition, X(t=0) = [S,E,I,R].
%     T - scalar
%       Determines the maximum extent of time the simulation computes,
%       modelling from t=0 to t=T.
%
%   Output Arguments
%     sample_success - (n x 1) matrix, elements ∈[0,1] 
%       The proportion of successful trials for a given paramater sample.


%% Test samples
testtimes = [0, 10, T];
n = size(samples,1);
sample_success = zeros(n,1);

for sample = 1:n
    % Establishing parameters to test
    alpha = samples(sample,1);
    beta = samples(sample,2);
    rho = samples(sample,3);
    for trial = 1:m 
        % Running simulation----------------------------------------------]
        [t, X] = gillespieSSA(alpha, beta, rho, X0, T); 
        S = X(1,:); % Extracting susceptible subpopulation

        % Extracting values for conditions--------------------------------]
        s = zeros(1,3); % Initialising s at t = testtimes
        if t(end) < testtimes(2),       s(2:3) = S(end); index = 1;
        elseif t(end) < testtimes(3),     s(3) = S(end); index = 1:2;
        else, index = 1:3; 
        end, s(index) = interp1(t', S', testtimes(index));

        % Testing parameter performance against given conditions----------]
        cond1 = s(1) - 20 <= s(3)   &&   s(3) <= s(2) - 10;
        cond2 = s(3) < s(1) - 100;
        if cond1 || cond2 % Count successful trials
            sample_success(sample) = sample_success(sample) + 1;
        end
    end % Finish trials for the sample
end % Finish testing all samples

sample_success = sample_success / m; % Proportion of successful trials

end