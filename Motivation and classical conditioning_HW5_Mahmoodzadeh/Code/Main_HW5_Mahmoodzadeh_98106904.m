%%
% Ali Akbar Mahmoodzadeh 98106904 HW5
%% Part 1: 
% Extiction:
clear; clc; close all;
num_trials = 200;
r = [ones(1, num_trials/2) zeros(1, num_trials/2)];
lr = 0.05;
w(1) = 0;
for itrial = 1: num_trials
    w(itrial + 1) = w(itrial) + lr * (r(itrial) - w(itrial));
end

figure;

plot(w, 'linewidth', 2);
title(['Extinction, Learning rate: ', num2str(lr)], 'interpreter', 'latex');
xlabel(['Trials'], 'interpreter', 'latex');
ylabel(['w'], 'interpreter', 'latex');
xlim([1, num_trials]);
grid on;

clear w1 w2 u1 u2 r



% === Partial:
num_trials = 6000;
alpha_vec = [0.2  0.6  0.9];

figure;
for ialpha = 1: length(alpha_vec)
    r = zeros(1, num_trials);
    alpha = alpha_vec(ialpha);
    randVec = randperm(num_trials);
    r(randVec(1: num_trials *alpha)) = 1;
    
    lr = 0.05;
    w(1) = 0;
    for itrial = 1: num_trials
        w(itrial + 1) = w(itrial) + lr * (r(itrial) - w(itrial));
    end
    disp(w(end))
    plot(w, 'linewidth', 1);
    grid on;
    hold all
end
title(['Partial, Learning rate: ', num2str(lr)], 'interpreter', 'latex');
xlabel(['Trials'], 'interpreter', 'latex');
ylabel(['w'], 'interpreter', 'latex');
xlim([1, num_trials]);
yline(0.2, 'k', 'linewidth', 3);

yline(0.6, 'k', 'linewidth', 3);

yline(0.9, 'k', 'linewidth', 3);
legend('alpha = 0.2', 'alpha = 0.5', 'alpha = 0.6', 'alpha = 0.8', 'alpha = 0.9', 'interpreter', 'latex');

clear w1 w2 u1 u2 r



% blocking 
num_trials = 400;
r = [ones(1, num_trials/2)];
lr = 0.05;
w1(1) = 0;
% pre-trainning
for itrial = 1: num_trials/2
    w1(itrial + 1) = w1(itrial) + lr * (r(itrial) - w1(itrial));
end
% training
w1new(1) = w1(end);
w2(1) = 0;
for itrial = 1: num_trials/2
    w1new(itrial + 1) = w1new(itrial) + lr * (r(itrial) - w1new(itrial) - w2(itrial));
    w2(itrial + 1) = w2(itrial) + lr * (r(itrial) - w1new(itrial) - w2(itrial));
end
figure;
plot(w1, 'linewidth', 2);
xlim([1, num_trials/2]);
grid on;
hold all;
plot(w2, 'linewidth', 2);
xlim([1, num_trials/2]);
grid on;
title(['Blocking, Learning rate: ', num2str(lr)], 'interpreter', 'latex');
xlabel(['Trials'], 'interpreter', 'latex');
ylabel(['w'], 'interpreter', 'latex');
legend('w1', 'w2',  'interpreter', 'latex');

clear w1 w2 u1 u2 r

% Inhibitory:
num_trials = 300;
lr = 0.05;

alpha = 0.2;
r = zeros(1, num_trials);
randVec = randperm(num_trials);
r(randVec(1: num_trials * alpha)) = 1;

u1 = ones(1, num_trials);
u2 = zeros(1, num_trials);
u2(randVec(num_trials * alpha + 1: end)) = 1;

w1(1) = 0;
w2(1) = 0;
for itrial = 1: num_trials
    if (r(itrial) == 0)
        w1(itrial + 1) = w1(itrial) + lr * (r(itrial) - w1(itrial) * u1(itrial)...
            - w2(itrial) * u2(itrial));
        w2(itrial + 1) = w2(itrial) + lr * (r(itrial) - w1(itrial) * u1(itrial)...
            - w2(itrial) * u2(itrial));
    else
        w1(itrial + 1) = w1(itrial) + lr * (r(itrial) - w1(itrial) * u1(itrial));
        w2(itrial + 1) = w2(itrial);
    end   
end
disp(w1(end))
disp(w2(end))

figure;
plot(w1, 'linewidth', 2);
xlim([1, num_trials]);
hold all;
plot(w2, 'linewidth', 2);
xlim([1, num_trials]);
grid on;
title(['Inhibitory (Learning rate: ', num2str(lr)],...
    'interpreter', 'latex');
xlabel(['Trials'], 'interpreter', 'latex');
ylabel(['w'], 'interpreter', 'latex');
legend('w1', 'w2',  'interpreter', 'latex');

clear w1 w2 u1 u2 r

% overshadow:
num_trials = 100;
lr1 = 0.05;
lr2 = 0.05;
r = ones(1, num_trials);
w1(1) = 0;
w2(1) = 0;
for itrial = 1: num_trials
    w1(itrial + 1) = w1(itrial) + lr1 * (r(itrial) - w1(itrial) - w2(itrial));
    w2(itrial + 1) = w2(itrial) + lr2 * (r(itrial) - w1(itrial) - w2(itrial));
end
disp(w1(end));
disp(w2(end));

figure;
plot(w1, 'linewidth', 2);
xlim([1, num_trials]);
grid on;
hold all;
plot(w2,'--r', 'linewidth', 2);
xlim([1, num_trials]);
grid on;
title(['Overshadow (Learning rate 1: ', num2str(lr1), ', Learning rate 2: ', ...
    num2str(lr2),')' ], 'interpreter', 'latex');
xlabel(['Trials'], 'interpreter', 'latex');
ylabel(['w'], 'interpreter', 'latex');
legend('w1', 'w2',  'interpreter', 'latex');

clear w1 w2 u1 u2 r


num_trials = 50;
lr1 = 0.05;
lr2 = 0.2;
r = ones(1, num_trials);
w1(1) = 0;
w2(1) = 0;
for itrial = 1: num_trials
    w1(itrial + 1) = w1(itrial) + lr1 * (r(itrial) - w1(itrial) - w2(itrial));
    w2(itrial + 1) = w2(itrial) + lr2 * (r(itrial) - w1(itrial) - w2(itrial));
end
disp(w1(end));
disp(w2(end));
figure;
plot(w1, 'linewidth', 2);
xlim([1, num_trials]);
grid on;
hold all;
plot(w2,'--r', 'linewidth', 2);
xlim([1, num_trials]);
grid on;
title(['Overshadow (Learning rate 1: ', num2str(lr1), ', Learning rate 2: ', ...
    num2str(lr2),')' ], 'interpreter', 'latex');
xlabel(['Trials'], 'interpreter', 'latex');
ylabel(['w'], 'interpreter', 'latex');
legend('w1', 'w2',  'interpreter', 'latex');

clear w1 w2 u1 u2 r

%% PART 2: KALMAN FILTER

clear all; close all; clc;

nTrials = 22;


s = rng;
v1 = normrnd(0,0.1,1,nTrials);
v2 = normrnd(0,0.1,1,nTrials);
w1 = zeros(1,nTrials);
w1(1) = 1;
w2 = zeros(1,nTrials);
w2(1) = 1;
for i = 2:nTrials
    w1(i) = w1(i-1) + v1(i-1);
    w2(i) = w2(i-1) + v2(i-1);
end



% Blocking
r = ones(1,nTrials);
w0 = [0;0];
C = [ones(1,nTrials);zeros(1,nTrials/2),ones(1,nTrials/2)];
W_noise = eye(2)*0.01; 
tau = 0.7; 
sigma0 = eye(2)*0.6;
[w,sigma] = KALMAN(W_noise,w0,nTrials,tau,sigma0,C,r);
sigma1 = [];
sigma2 = [];
for i = 1:nTrials
    tmp = sigma{i};
    sigma1 = [sigma1, tmp(1,1)];
    sigma2 = [sigma2, tmp(2,2)];
end

figure;
subplot(2,1,1)
plot(0:nTrials-1,w(1,:),'k','LineWidth',2)
hold on
plot(0:nTrials-1,w(2,:),'--r','LineWidth',2)
hold on
xline(10,':k');
xlim([0 20])
ylim([0 1.2])
title("Blocking      mean",'Interpreter','latex')
ylabel("W",'Interpreter','latex')
legend("W_1","W_2",'Interpreter','latex')
grid on;

subplot(2,1,2)
plot(0:nTrials-1,sigma1,'k','LineWidth',2)
hold on 
plot(10:nTrials-1,sigma2(11:end),'--r','LineWidth',2)
hold on
xline(10,':k');
xlim([0 20])
ylim([0 1])
title("Blocking      Variance",'Interpreter','latex')
ylabel("$\sigma^2(t)$",'Interpreter','latex')
legend("${\sigma{^2}}_1$","${\sigma{^2}}_2$",'Interpreter','latex')
grid on;

% Unblocking
r = [ones(1,nTrials/2),2*ones(1,nTrials/2)];
w0 = [0;0];
C = [ones(1,nTrials);zeros(1,nTrials/2),ones(1,nTrials/2)];
W_noise = eye(2)*0.01;
tau = 0.6;
sigma0 = eye(2)*0.6;
[w,sigma] = KALMAN(W_noise,w0,nTrials,tau,sigma0,C,r);
sigma1 = [];
sigma2 = [];
for i = 1:nTrials
    tmp = sigma{i};
    sigma1 = [sigma1, tmp(1,1)];
    sigma2 = [sigma2, tmp(2,2)];
end

figure;
subplot(2,1,1)
plot(0:nTrials-1,w(1,:),'k','LineWidth',2)
hold on
plot(0:nTrials-1,w(2,:),'--r','LineWidth',2)
hold on
xline(10,':k');
xlim([0 20])
ylim([0 1.2])
title("Unblocking      mean",'Interpreter','latex')
ylabel("W",'Interpreter','latex')
legend("W_1","W_2",'Interpreter','latex')
grid on;

subplot(2,1,2)
plot(0:nTrials-1,sigma1,'k','LineWidth',2)
hold on 
plot(10:nTrials-1,sigma2(11:end),'--r','LineWidth',2)
hold on
xline(10,':k');
xlim([0 20])
ylim([0 1])
title("Unblocking      Variance",'Interpreter','latex')
ylabel("$\sigma^2(t)$",'Interpreter','latex')
legend("${\sigma{^2}}_1$","${\sigma{^2}}_2$",'Interpreter','latex')
grid on;


% Backward Blocking
r = ones(1,nTrials);
w0 = [0;0];
C = [ones(1,nTrials);ones(1,nTrials/2),zeros(1,nTrials/2)];
W_noise = eye(2)*0.02;
tau = 1.2;
sigma0 = eye(2)*0.6;
[w,sigma] = KALMAN(W_noise,w0,nTrials,tau,sigma0,C,r);
sigma1 = [];
sigma2 = [];
for i = 1:nTrials
    tmp = sigma{i};
    sigma1 = [sigma1, tmp(1,1)];
    sigma2 = [sigma2, tmp(2,2)];
end

figure;
subplot(2,1,1)
plot(0:nTrials-1,w(1,:),'k','LineWidth',2)
hold on
plot(0:nTrials-1,w(2,:),'--r','LineWidth',2)
hold on
xline(10,':k');
xlim([0 20])
ylim([0 1.2])
title("Backward Blocking      Mean",'Interpreter','latex')
ylabel("W",'Interpreter','latex')
legend("W_1","W_2",'Interpreter','latex')
grid on;

subplot(2,1,2)
plot(0:nTrials-1,sigma1,'k','LineWidth',2)
hold on 
plot(10:nTrials-1,sigma2(11:end),'--r','LineWidth',2)
hold on
xline(10,':k');
xlim([0 20])
ylim([0 1])
title("Backward Blocking      Variance",'Interpreter','latex')
ylabel("$\sigma^2(t)$",'Interpreter','latex')
legend("${\sigma{^2}}_1$","${\sigma{^2}}_2$",'Interpreter','latex')
grid on;



%% Figure 1B
clear; clc; close all;
% rng('default');
num_trials = 20;
w1(1) = 1;
w2(1) = 1;
for itrial = 1: num_trials - 1
    w1(itrial + 1) =  w1(itrial) + normrnd(0, 0.1);
    w2(itrial + 1) =  w2(itrial) + normrnd(0, 0.1);
end
figure;
plot(w1, 'black');
hold all
plot(w2, '--b', 'color', 'black');
xlabel(['Trials'], 'interpreter', 'latex');
ylabel(['W'], 'interpreter', 'latex');
legend('W1', 'W2', 'interpreter', 'latex');
title('Drift', 'interpreter', 'latex');
grid on;








%% 5:
clear; clc; close all;

num_trials = 20;
r = [ones(1, num_trials/2), -1 * ones(1, num_trials/2)];
u = [ones(1, num_trials)]';
w0 = 0;
sigma0 = 0.6;
W = 0.01;
tau = 0.47;
[w, sigma] = myKalman(r, u, w0, sigma0, W, tau, num_trials, 0);
sigmavec = [];
for i = 1: num_trials
    sigmavec = [sigmavec, sigma{i}];
end
figure;
subplot(2, 1, 1)
plot(w ,'k','LineWidth',2);
ylabel('w', 'interpreter', 'latex');
xlabel('Trials', 'interpreter', 'latex');
grid on;
subplot(2, 1, 2)
plot(sigmavec ,'k','LineWidth',2)
ylabel('$\sigma^2$', 'interpreter', 'latex');
xlabel('Trials', 'interpreter', 'latex');
grid on;
sgtitle('s1 $\rightarrow$ r, s1 $\rightarrow$ -r subject to: (W = 0.01, T = 0.4, $\Sigma_0$ = 0.6I)', ...
    'interpreter', 'latex');











%% FUNCTIONS

function [w, sigma] = myKalman(r, u, w0, sigma0, W, tau, num_trials, plot_flag)

    w(:, 1) = w0;
    sigma{1} = sigma0;
    for itrial = 1: num_trials - 1
        if (itrial == num_trials/2 + 1 && size(u,2) > 1)
            sigma{itrial}(size(u,2)^2) = sigma0(size(u,2)^2);
        end
        sigma_prd{itrial} = sigma{itrial} + W;
        G = sigma_prd{itrial} * u(itrial, :)' * inv(u(itrial, :) * sigma_prd{itrial}...
            * u(itrial, :)' + tau^2);
        sigma{itrial + 1} = sigma_prd{itrial} - G * u(itrial, :) * sigma_prd{itrial};
        w(:, itrial + 1) = w(:, itrial) + G * (r(itrial) - u(itrial, :) * w(:, itrial));
    end
    
    if plot_flag
    
        figure;
        subplot(2, 1, 1)
        plot([0: num_trials-1], w(1, :), 'color', 'black');
        hold all;

        plot([0: num_trials-1], w(2, :), '--', 'color', 'black');
                
        xline(num_trials/2, 'LineStyle', ':', 'color', 'black');
        xlabel('Trials', 'interpreter', 'latex');
        ylabel('w', 'interpreter', 'latex');
        legend('$w_1$', '$w_2$',  'interpreter', 'latex');
        
        sigma1 = [];
        sigma2 = [];
        for itrial = 1: num_trials
            sigma1 = [sigma1, sigma{itrial}(1)];
            sigma2 = [sigma2, sigma{itrial}(4)];
        end
        subplot(2, 1, 2)
        plot([0: num_trials-1], sigma1, 'color', 'black');
        hold all;
        plot([num_trials/2: num_trials-1], sigma2(num_trials/2 + 1: end), '--', 'color', 'black');

        xline(num_trials/2, 'LineStyle', ':', 'color', 'black');
        xlabel('Trials', 'interpreter', 'latex');
        ylabel('$\sigma^2$', 'interpreter', 'latex');
        legend('$\sigma_1$', '$\sigma_2$',  'interpreter', 'latex');
    
    end
end











