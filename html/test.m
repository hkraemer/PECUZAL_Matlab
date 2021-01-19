clear
set(0, 'DefaultAxesFontSize', 14, 'DefaultLineLineWidth', 1) % increase font size and line width for the figures
  
% set parameters for standard Roessler
a = 0.2;
b = 0.2;
c = 5.7;

% Set time series length and sampling
N = 5000;
dt = 0.2;
transient_time = 2500;
tf = (N+transient_time) * dt; 
t = 0:dt:tf;

roe = @(t,x) [-x(2)-x(3); x(1)+a*x(2); b+x(3)*(x(1)-c)];

% solve the ODE
x0 = [1.; .5; 0.5]; % initial conditions  
[t, x] = ode45(roe, t, x0);

% remove transients
x(1:transient_time+1,:) = [];
t = t(1:end-transient_time-1);

%%
y = x(:,2);  % bind y-component
mi = mutualinformation(y,50); % compute MI

figure('Units', 'normalized', 'Position', [.3 .3 .7 .5])

subplot(1,2,1)
plot(t(1:1000), y(1:1000))
title('y-component of Roessler time series')
xlabel('time')
grid on

subplot(1,2,2)
plot(mi(:,1), mi(:,2), '.-', 'MarkerSize',15)
xlabel('time delay \tau')
ylabel('mutual information [nats]')
title('Mutual information for y-component of Roessler time series')
grid on

%%
tic
[Y_reconstruct, tau_vals, ts_vals, Ls, E] = pecuzal_embedding(y, 0:100, 'theiler', 7, 'econ', true);
toc