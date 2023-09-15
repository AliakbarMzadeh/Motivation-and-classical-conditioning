%% Fig 3


Gamma = 3;

figure;
set(gcf,'Color',[1 1 1]);

subplot(3,1,1);
w0 = [0:-0.02:-1 -3:-0.02:-3.5];
r = w0' + normrnd(0,0.2,length(w0),1);
w0 = w0' + normrnd(0,0.05,length(w0),1);
plot(w0,'k','LineWidth',2);
grid on;
hold on
scatter(1:length(w0),r,'x','g','LineWidth',1);
ylim([-4 4])
grid on;
legend({'r(t)' 'w(t)'}, 'interpreter', 'latex');
box off
subplot(3,1,2);
plot(w0,'k','LineWidth',2);
grid on;
hold on
scatter(1:length(w0),r,'x','b','LineWidth',1);
grid on;
ylim([-4 4])


M = 0.49; 
P = 0.01; 
w = 0;
sig = 0.6;
u = 1;
w_list = [];
sig_list = [];
beta_list = [];
for i=1:length(r)
    y = r(i);
    Beta = (y - u*w)^2/(u'*sig*u+M);
    if(Beta > Gamma)
        sig = 1;
    end
    w_n = w + sig * u' ./ (u*sig*u'+M) * (y - u*w);
    sig_n = sig - sig * u' ./ (u*sig*u'+M) * u * sig + eye(1) * P;
    w = w_n;
    sig = sig_n;
    beta_list = [beta_list Beta];
    w_list = [w_list w];
    sig_list = [sig_list sig];
end
scatter(1:length(w0),w_list,'o','r','LineWidth',1);
grid on;
legend({'r(t)' 'w(t)' 'w^(t)'}, 'interpreter', 'latex');
box off
subplot(3,1,3);
plot(beta_list,'-')
grid on;
hold on
plot([0,length(w0)],[Gamma,Gamma])
grid on;
legend({'NE' 'Gamma'}, 'interpreter', 'latex');
box off

