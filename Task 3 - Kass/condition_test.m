function sample_success = condition_test(samples, m, X0, T)
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