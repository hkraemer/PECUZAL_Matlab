%% Test PECUZAL functionality

clear

% Test case for univariate example
addpath('data')
addpath('../.')
data = load('lorenz_pecora_uni_x.csv');
data = data(1:500);
theiler = 21;
Tmax = 100;

taus = 0:Tmax;

[~, tau_vals, ts_vals, Ls, ~] = pecuzal_embedding(data, taus, 'theiler', theiler);

assert(-0.74 < sum(Ls))
assert(sum(Ls) < -0.73)
assert(tau_vals(2) == 21)
assert(tau_vals(3) == 78)
assert(tau_vals(4) == 6)
assert(length(ts_vals) == 4)

[~, tau_vals, ~, ~, ~] = pecuzal_embedding(data, taus, 'theiler', theiler, 'L_thres', 0.2);

assert(length(tau_vals) == 3)

[~, tau_vals, ts_vals, Ls, ~] = pecuzal_embedding(data, taus, 'theiler', theiler, 'econ', true);

assert(-0.733 < sum(Ls))
assert(sum(Ls) < -0.73)
assert(tau_vals(2) == 21)
assert(tau_vals(3) == 78)
assert(tau_vals(4) == 6)
assert(length(ts_vals) == 4)

rng(1)
[~, tau_vals, ts_vals, Ls, ~] = pecuzal_embedding(data, taus, 'theiler', theiler, 'econ', true, 'sample_size', 0.8);

assert(-0.691 < sum(Ls))
assert(sum(Ls) < -0.692)
assert(tau_vals(2) == 20)
assert(tau_vals(3) == 76)
assert(tau_vals(3) == 10)
assert(length(ts_vals) == 4)


% Test case for multivariate example
cd('./tests/data')
data = load('lorenz_pecora_multi.csv');
cd('../..')
data1 = data(1:500,1:2);
data2 = data(1:7000,1:2);
data3 = data(1:8000,1:2);
data4 = data(1:9000,1:2);
data5 = data(1:10000,1:2);
theiler = 15;
Tmax = 100;
taus = 0:Tmax;
Tw = 10;

% uzal cost
L1 = uzal_cost(data2, theiler, 3, Tw, 1);
L2 = uzal_cost(data3, theiler, 3, Tw, 1);
L3 = uzal_cost(data4, theiler, 3, Tw, 1);
L4 = uzal_cost(data5, theiler, 3, Tw, 1);

delta_L2 = uzal_cost_pecuzal(data3(:,1), data3, Tw, 'theiler', theiler);
delta_L3 = uzal_cost_pecuzal(data4(:,1), data4, Tw, 'theiler', theiler);
delta_L4 = uzal_cost_pecuzal(data5(:,1), data5, Tw, 'theiler', theiler);

assert(L1 < L2)
assert(L2 < L3)
assert(L4 > L1)
assert(L3 < -0.475)
assert(L3 > -0.485)
assert(delta_L2 < -0.766)
assert(delta_L3 < -0.766)
assert(delta_L4 < -0.766)
assert(delta_L2 > -0.771)
assert(delta_L3 > -0.771)
assert(delta_L4 > -0.771)

samp = 0.5;
rng(1)
delta_L22 = uzal_cost_pecuzal(data3(:,1), data3, Tw, 'theiler', theiler, 'samplesize', samp);

assert(delta_L22 < -0.752)
assert(delta_L22 > -0.753)

% pecuzal
[Y, tau_vals, ts_vals, Ls, ~] = pecuzal_embedding(data1, taus, 'theiler', theiler);

assert(size(Y,2)==3)
assert(tau_vals(1) == 0)
assert(tau_vals(2) == 25)
assert(tau_vals(3) == 81)
assert(ts_vals(1) == 2)
assert(ts_vals(2) == 1)
assert(ts_vals(3) == 1)

assert(sum(Ls) < -0.6211)

[Y, tau_vals, ts_vals, Ls, ~] = pecuzal_embedding(data1, taus, 'theiler', theiler, 'econ', true);

assert(size(Y,2)==3)
assert(tau_vals(1) == 0)
assert(tau_vals(2) == 25)
assert(tau_vals(3) == 81)
assert(ts_vals(1) == 2)
assert(ts_vals(2) == 1)
assert(ts_vals(3) == 1)

assert(sum(Ls) < -0.62)

rng(1)
[Y, tau_vals, ts_vals, Ls, ~] = pecuzal_embedding(data1, taus, 'theiler', theiler, 'econ', true, 'sample_size', 0.5);

assert(size(Y,2)==5)
assert(tau_vals(1) == 0)
assert(tau_vals(2) == 0)
assert(tau_vals(3) == 78)
assert(tau_vals(4) == 37)
assert(tau_vals(5) == 24)
assert(ts_vals(1) == 2)
assert(ts_vals(2) == 1)
assert(ts_vals(3) == 1)
assert(ts_vals(4) == 1)
assert(ts_vals(5) == 1)

assert(sum(Ls) < -0.759)
