function prob_table = get_binomial_table(p,alpha,delta_max)
% for the Pecora-FNN method we need to estimate how many points from the
% delta neighborhood need to fall outside the epsilon neighborhood to
% reject the Null Hypothesis at a certain significance level alpha. Here we
% compute a table storing these numbers of points for a binomial parameter
% `p`, a significance level `alpha` and a range of trials (`n`'s) from 
% 8:`delta_max`.

% K.H.Kraemer, May 2020

% set range for trials
n_min = 8;
n_max = delta_max;

prob_table = zeros(length(n_min:n_max),3);
prob_table(:,1) = n_min:n_max;

% loop over the different n values
cnt = 1;
for n = n_min:n_max
    prob_table(cnt,2) = binoinv(1-alpha,n,p);
    cnt = cnt + 1;
end
prob_table(:,3)= prob_table(:,1)-prob_table(:,2);

end