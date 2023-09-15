
%% Functions 


function [w,sigma] = KALMAN(W_noise,w0,nTrials,tau,sigma0,C,r)
    w = zeros(size(C,1),nTrials);
    sigma = cell(1,nTrials);
    w(:,1) = w0;
    sigma{1} = sigma0;
    if size(C,1) > 1
        ind = find(C(2,:));
        ind = ind(1)-1;
    end
    for i = 2:nTrials
        % Prediction
        wp = w(:,i-1);
        sigmap = sigma{i-1} + W_noise;
        
        % Update
        G = sigmap * C(:,i) * (C(:,i)' * sigmap * C(:,i) + tau^2)^-1;
        sigma{i} = sigmap - G*C(:,i)'*sigmap;
        (r(i) - C(:,i)'*w(:,i-1));
        w(:,i) = w(:,i-1) + G * (r(i) - C(:,i)'*w(:,i-1));
        
        if size(C,1) > 1
            if (i == ind)
                tmp = sigma{i};
                tmp(2,2) = 0.6;
                sigma{i} = tmp;
            end
        end
    end
end










