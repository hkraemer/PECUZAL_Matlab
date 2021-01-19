clear, clc

data = load('data/lorenz_pecora_uni_x.csv');
data = data(1:500);
theiler = 21;
Tmax = 100;

%%
s = round(data,4);

%%
taus = [0,58,12];
js = [1,1,1];
x = s;
for m = 1:length(taus)
    if m == 1
        Y_old = embed(x(:,js(m)), m, taus(m));
    else
        Y_old = embed2(Y_old, x(:,js(m)), taus(m));
    end
end

fiducial_point = 1;
theiler = 21;
[distances, ~] = all_distances(Y_old(fiducial_point,:),...
                              Y_old(1:end-100,:), 'euc');
neighbours = 13;

[~,ind] = sort(distances);

% temporal neighbours within Theiler window
idx = max(1, fiducial_point - theiler):min(length(Y_old)-100, fiducial_point + theiler);
% remove these neighbours from index list for distances
ind(ismember(ind, idx)) = [];
% construct neighbourhood
NN_idxs = ind(1:neighbours)';
NN_dist = distances(ind(1:neighbours));


%%
% eps = pecora_embedding_cycle(data, [0], [1], 0:100, 1.0, 21);
eps = pecora_embedding_cycle(s, [0,58,12], [1,1,1], 0:100, 1.0, 21);


figure
plot(eps)
grid on

%%
taus = [0,2];
js = [1,1];
x = s;
for m = 1:length(taus)
    if m == 1
        Y_old = embed(x(:,js(m)), m, taus(m));
    else
        Y_old = embed2(Y_old, x(:,js(m)), taus(m));
    end
end

%%
fiducial_point = 1;
[distances, ~] = all_distances(Y_old(fiducial_point,:),...
                              Y_old(1:end-100,:), 'euc');
neighbours = 13;

[~,ind] = sort(distances);

% temporal neighbours within Theiler window
idx = max(1, fiducial_point - theiler):min(398, fiducial_point + theiler);
% remove these neighbours from index list for distances
ind(ismember(ind, idx)) = [];
% construct neighbourhood
NN_idxs = ind(1:neighbours)';
NN_dist = distances(ind(1:neighbours));
%%
