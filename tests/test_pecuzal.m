%% Test PECUZAL functionality

clear

% Test case for univariate example
data = load('./data/lorenz_pecora_uni_x.csv');
data = data(1:500);
theiler = 21;
Tmax = 100;

taus = 0:Tmax;

%%

[~, tau_vals, ts_vals, Ls, ~] = pecuzal_embedding(data, taus, 'theiler', theiler);

assert(-0.6134 < sum(Ls))
assert(sum(Ls) < -0.6133)
assert(tau_vals(2) == 21)
assert(tau_vals(3) == 13)
assert(tau_vals(4) == 78)

assert(length(ts_vals) == 4)

[~, tau_vals, ~, ~, ~] = pecuzal_embedding(data, taus, 'theiler', theiler, 'L_thres', 0.2);

assert(length(tau_vals) == 2)


[~, tau_vals, ts_vals, Ls, ~] = pecuzal_embedding(data, taus, 'theiler', theiler, 'econ', true);

assert(-0.6078 < sum(Ls))
assert(sum(Ls) < -0.6077)
assert(tau_vals(2) == 21)
assert(tau_vals(3) == 13)
assert(tau_vals(4) == 78)

assert(length(ts_vals) == 4)

%% Test case for multivariate example
data = load('data/lorenz_pecora_multi.csv');
data1 = data(1:500,1:2);
data2 = data(1:1000,1:2);
data3 = data(1:2000,1:2);
theiler = 15;
Tmax = 100;
taus = 0:Tmax;
Tw = 10;

L1 = uzal_cost(data1, theiler, 3, Tw, 1);
L2 = uzal_cost(data2, theiler, 3, Tw, 1);
L3 = uzal_cost(data3, theiler, 3, Tw, 1);

delta_L1 = uzal_cost_pecuzal(data1(:,1), data1, Tw, 'theiler', theiler);
delta_L2 = uzal_cost_pecuzal(data2(:,1), data2, Tw, 'theiler', theiler);
delta_L3 = uzal_cost_pecuzal(data3(:,1), data3, Tw, 'theiler', theiler);

%%
[Y, tau_vals, ts_vals, Ls, ~] = pecuzal_embedding(data1, taus, 'theiler', theiler);

assert(size(Y,2)==2)
assert(tau_vals(1) == 0)
assert(tau_vals(2) == 0)
assert(ts_vals(1) == 2)
assert(ts_vals(2) == 1)

assert(sum(Ls) < -0.5505736)

[Y, tau_vals, ts_vals, Ls, ~] = pecuzal_embedding(data1, taus, 'theiler', theiler, 'econ', true);

assert(size(Y,2)==2)
assert(tau_vals(1) == 0)
assert(tau_vals(2) == 0)
assert(ts_vals(1) == 2)
assert(ts_vals(2) == 1)

assert(sum(Ls) < -0.544942)


