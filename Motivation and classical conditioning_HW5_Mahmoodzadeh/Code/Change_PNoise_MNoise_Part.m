
%% Q2 +++++++======
% Fig 1 Blocking


% Define parameter combinations
param_combinations = [
    0.01 0.5;
    0.1 0.5;
    0.5 0.5;
    0.8 0.5;
    0.5 0.01;
    0.5 0.1;
    0.5 0.8;
    0.5 0.9
];

n_Trials = 10;
num_trials = 20;

for i = 1:size(param_combinations, 1)
    
    M = param_combinations(i, 1);
    P = param_combinations(i, 2);

    w = [0; 0];
    sig = [0.6 0; 0 0.6];
    u = [1 0];
    r = 1;
    w_list = w;
    sig_list = [sig(1,1);sig(2,2)];

    for j = 1:n_Trials
        y = r + M;
        w_n = w + sig * u' ./ (u*sig*u'+M) * (y - u*w);
        sig_n = sig - sig * u' ./ (u*sig*u'+M) * u * sig + eye(2) * P;
        w = w_n;
        sig = sig_n;
        w_list = [w_list w];
        sig_list = [sig_list [sig(1,1);sig(2,2)]];
    end

    u = [1 1];
    for j = 1:n_Trials
        y = r + M;
        w_n = w + sig * u' ./ (u*sig*u'+M) * (y - u*w);
        sig_n = sig - sig * u' ./ (u*sig*u'+M) * u * sig + eye(2) * P;
        w = w_n;
        sig = sig_n;
        w_list = [w_list w];
        sig_list = [sig_list [sig(1,1);sig(2,2)]];
    end

    % Plotting
    figure;
    set(gcf,'Color',[1 1 1]);

    subplot(2,1,1);
    plot(0:(2*n_Trials),w_list(1,:),'k','LineWidth',2);
    hold on
    plot(0:(2*n_Trials),w_list(2,:),'--r','LineWidth',2);
xline(num_trials/2,'k','LineWidth',2);
    xlabel('Trial', 'interpreter', 'latex');
    ylabel('W', 'interpreter', 'latex');
    legend({'W1', 'W2'}, 'interpreter', 'latex');
    grid on;
    box off
    

    title(['Measurement Noise: ' num2str(M) ', Process Noise: ' num2str(P)], 'interpreter', 'latex');

    subplot(2,1,2);
    plot(0:(2*n_Trials),sig_list(1,:),'k','LineWidth',2);
    hold on
    plot(n_Trials:(2*n_Trials),sig_list(2,n_Trials+1:end),'--r','LineWidth',2);
    box off
   xline(num_trials/2,'k','LineWidth',2);
    xlabel('Trial', 'interpreter', 'latex');
    ylabel("$\sigma^2(t)$",'Interpreter','latex')
        
    
    legend("${\sigma{^2}}_1$","${\sigma{^2}}_2$",'Interpreter','latex')
    grid on;

   
end

%% Q2++++++++==========
% Fig 1 Unblocking


% Define parameter combinations
param_combinations = [
    0.0 0.5;
    0.1 0.5;
    0.5 0.5;
    0.8 0.5;
    0.5 0.0;
    0.5 0.1;
    0.5 0.8;
    0.5 0.9
];
num_trials = 20;

n_Trials = 10;

for i = 1:size(param_combinations, 1)
    
    P = param_combinations(i, 1);
    M = param_combinations(i, 2);

    w = [0; 0];
    sig = [0.6 0; 0 0.6];
    u = [1 0];
    r = 1;
    w_list = w;
    sig_list = [sig(1,1);sig(2,2)];

    for j = 1:n_Trials
        y = r + M;
        w_n = w + sig * u' ./ (u*sig*u'+M) * (y - u*w);
        sig_n = sig - sig * u' ./ (u*sig*u'+M) * u * sig + eye(2) * P;
        w = w_n;
        sig = sig_n;
        w_list = [w_list w];
        sig_list = [sig_list [sig(1,1);sig(2,2)]];
    end

    u = [1 1];
    r = 2;
    for j = 1:n_Trials
        y = r + M;
        w_n = w + sig * u' ./ (u*sig*u'+M) * (y - u*w);
        sig_n = sig - sig * u' ./ (u*sig*u'+M) * u * sig + eye(2) * P;
        w = w_n;
        sig = sig_n;
        w_list = [w_list w];
        sig_list = [sig_list [sig(1,1);sig(2,2)]];
    end

    % Plotting
    figure;
    set(gcf,'Color',[1 1 1]);

    subplot(2,1,1);
    plot(0:(2*n_Trials),w_list(1,:),'k','LineWidth',2);
    hold on
    plot(0:(2*n_Trials),w_list(2,:),'--r','LineWidth',2);
    xline(num_trials/2,'k','LineWidth',2);
    xlabel('Trial', 'interpreter', 'latex');
    ylabel('W', 'interpreter', 'latex');
    legend({'W1', 'W2'}, 'interpreter', 'latex');
    grid on;
    box off
    title(['Measurement Noise: ' num2str(M) ', Process Noise: ' num2str(P)],'Interpreter','latex');

    subplot(2,1,2);
    plot(0:(2*n_Trials),sig_list(1,:),'k','LineWidth',2);
    hold on
    plot(n_Trials:(2*n_Trials),sig_list(2,n_Trials+1:end),'--r','LineWidth',2);
    box off
    xline(num_trials/2,'k','LineWidth',2);
    xlabel('Trial', 'interpreter', 'latex');
    ylabel("$\sigma^2(t)$",'Interpreter','latex')
        
    
    legend("${\sigma{^2}}_1$","${\sigma{^2}}_2$",'Interpreter','latex')
    grid on;

  
end







%% Q2++++++++++======
% Fig 1 Backwrd Blocking
P_Malues = [0, 0.1, 0.5, 0.8];
M_Malues = [0.5, 0.5, 0.5, 0.5];

for i = 1:numel(P_Malues)
    M = M_Malues(i);
    P = P_Malues(i);
    
    w = [0; 0];
    sig = [0.6 0; 0 0.6];
    n_Trials = 10;

    u = [1 1];
    r = 1;
    w_list = w;
    sig_list = [sig(1,1);sig(2,2);sig(1,2);sig(2,1)];
    for j = 1:n_Trials
        y = r + M;
        w_n = w + sig * u' ./ (u*sig*u'+M) * (y - u*w);
        sig_n = sig - sig * u' ./ (u*sig*u'+M) * u * sig + eye(2) * P;
        w = w_n;
        sig = sig_n;
        w_list = [w_list w];
        sig_list = [sig_list [sig(1,1);sig(2,2);sig(1,2);sig(2,1)]];
    end
    u = [1 0];
    for j = 1:n_Trials
        y = r + M;
        w_n = w + sig * u' ./ (u*sig*u'+M) * (y - u*w);
        sig_n = sig - sig * u' ./ (u*sig*u'+M) * u * sig + eye(2) * P;
        w = w_n;
        sig = sig_n;
        w_list = [w_list w];
        sig_list = [sig_list [sig(1,1);sig(2,2);sig(1,2);sig(2,1)]];
    end

    figure;
    set(gcf,'Color',[1 1 1]);

    subplot(1,3,1);
    w1 = -1:0.02:2;
    w2 = -1:0.02:2;
    [w1,w2] = meshgrid(w1,w2);
    X = [w1(:) w2(:)];
    y = mMnpdf(X,w_list(:,1)',[sig_list(1,1) sig_list(3,1); sig_list(4,1) sig_list(2,1)]);
    y = reshape(y,length(w2),length(w1));
    surf(w1,w2,y)
    caxis([min(y(:))-0.5*range(y(:)),max(y(:))])
    axis([-3 3 -3 3 0 0.4])
    xlim([-1 2])
    ylim([-1 2])
    xlabel('w1')
    ylabel('w2')
    zlabel('Probability Density')
    Miew([0 90]);
    shading interp
    colormap(gray)
    title(['t = 1, P = ' num2str(P) ', M = ' num2str(M)]);
    hold on
    scatter3(w_list(1,1),w_list(2,1),0.3,'*');
    axis square
    
    subplot(1,3,2);
    w1 = -1:0.02:2;
    w2 = -1:0.02:2
    subplot(1,3,2);
    [w1,w2] = meshgrid(w1,w2);
    X = [w1(:) w2(:)];
    y = mMnpdf(X,w_list(:,9)',[sig_list(1,9) sig_list(3,9); sig_list(4,9) sig_list(2,9)]);
    y = reshape(y,length(w2),length(w1));
    surf(w1,w2,y)
    caxis([min(y(:))-0.5*range(y(:)),max(y(:))])
    axis([-3 3 -3 3 0 0.4])
    xlim([-1 2])
    ylim([-1 2])
    xlabel('w1')
    ylabel('w2')
    zlabel('Probability Density')
    Miew([0 90]);
    shading interp
    colormap(gray)
    title(['t = 9, P = ' num2str(P) ', M = ' num2str(M)]);
    hold on
    scatter3(w_list(1,9),w_list(2,9),50,'*');
    axis square
    
    subplot(1,3,3);
    [w1,w2] = meshgrid(w1,w2);
    X = [w1(:) w2(:)];
    y = mMnpdf(X,w_list(:,19)',[sig_list(1,19) sig_list(3,19); sig_list(4,19) sig_list(2,19)]);
    y = reshape(y,length(w2),length(w1));
    surf(w1,w2,y)
    caxis([min(y(:))-0.5*range(y(:)),max(y(:))])
    axis([-3 3 -3 3 0 0.4])
    xlim([-1 2])
    ylim([-1 2])
    xlabel('w1')
    ylabel('w2')
    zlabel('Probability Density')
    Miew([0 90]);
    shading interp
    colormap(gray)
    title(['t = 19, P = ' num2str(P) ', M = ' num2str(M)]);
    hold on
    scatter3(w_list(1,19),w_list(2,19),50,'*');
    axis square
  
end


