%% Test PECUZAL functionality

% clear

assert(4==4)
% % Test case for univariate example
% data = load('./data/lorenz_pecora_uni_x.csv');
% data = data(1:500);
% theiler = 21;
% Tmax = 100;
% K = 3;
% KNN = 14;
% taus = 0:Tmax;
% 
% addpath('../src')
% 
% [~, tau_vals, ts_vals, Ls, ~] = pecuzal_embedding(data, taus, 'theiler', theiler, 'sample_size', 1, 'K', K, 'KNN', KNN);

% assert(-1.68 < Ls(1) && Ls(1) < -1.67)
% assert(-1.72 < Ls(2) && Ls(2) < -1.71)
% assert(-1.71 < Ls(3) && Ls(3) < -1.70)
% 
% assert(tau_vals(2) == 58)
% assert(tau_vals(3) == 12)
% 
% assert(length(ts_vals) == 3)
% 
% % Test case for multivariate example
% addpath('../tests')
% data = load('./data/lorenz_pecora_multi.csv');
% data = data(1:500,:);
% theiler = 15;
% Tmax = 100;
% taus = 0:Tmax;
% addpath('../src')
% [Y, tau_vals, ts_vals, Ls, eps] = pecuzal_embedding(data, taus, 'theiler', theiler);
% 
% assert(tau_vals(1) == tau_vals(2) && tau_vals(1) == 0)
% 
% assert(ts_vals(1) == 2)
% assert(ts_vals(2) == 3)
% assert(length(ts_vals) == 2)
% 
% assert(-1.69 < Ls(1) && Ls(1) < -1.68)
% assert(-1.68 < Ls(2) && Ls(2) < -1.67)