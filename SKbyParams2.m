%% SK model
%rng(4)
seed=8;
rng(seed)

%% Set dimension, trails, beta, external field strength
N = 10;
T = 1;
beta = 1;
%lambdas = [0:0.05:1];%0:0.05:1;%[0.01,0.1,0.25,0.5,0.75,0.9,0.99];
lambdas = 0.5;
h = 0;

predRhos = sign(2.*lambdas-1 + 4.*(1-lambdas)./N);

% rhos = zeros(numel(Ns),numel(betas));
% rhohats = zeros(numel(Ns),numel(betas));

% r2s = zeros(numel(Ns),numel(betas));
% gammas = zeros(numel(Ns),numel(betas));

rhos = zeros(T,numel(lambdas));
rhohats = zeros(T,numel(lambdas));
r2s = zeros(T,numel(lambdas));
gammas = zeros(T,numel(lambdas));

for k=1:numel(lambdas)
    lambda = lambdas(k);
    for t=1:T
        % Couplings
        J = normrnd(0,1/sqrt(N),N,N);
        
        % Energies
        E = SKEnergy(N,J,h);
        % Rates
        W = zeros(2^N);
        for l=1:2^N
            for m=1:2^N
                W(l,m) = exp(beta*(lambda*(E(l)+E(m))-E(m)));
            end
        end

        % Hypercube or complete graph
        A = double(hypercube(N));
        
        % Exit rates
        qout = sum(W.*A,2);
        W = W.*A - diag(qout);
        
        %% Solve for p, q, psi
        p = pFromQ(W);
        p = p';
        psi = (p.*qout)./sum(p.*qout);
        
        rho = corr(log(p),-log(qout));
        rhohat = corr(log(psi),-log(qout));
        r2 = var(log(psi))/var(log(qout));
        gamma = 1+rhohat*sqrt(r2);

        %% Calculate interesting quantities
        % rhos(i,j,k) = rho;
        % rhohats(i,j,k) = rhohat;
        % r2s(i,j,k) = r2;
        % gammas(i,j,k) = gamma;

        rhos(t,k) = rho
        rhohats(t,k) = rhohat;
        r2s(t,k) = r2;
        gammas(t,k) = gamma;
    end
    k
end

%scatter(sqrt(r2s),rhohats,'black','filled')

% mrho2 = mean(rhos);
% mrhohats2 = mean(rhohats);
% mr2s2 = mean(r2s);
% mgammas2 = mean(gammas);
% 
% figure
% plot(lambdas,mrho2,'LineWidth',2)
% hold on;
% plot(lambdas,predRhos,'LineWidth',2)
% xlabel('parameter of Glauber dynamics')
% ylabel('correlation')
% legend('actual','predicted')
% hold off
% ylim([-1,1])

% figure
% plot(lambdas,mrhohats2,'LineWidth',2)

% vrho2 = var(rhos,[],1);
% vrhohat = var(rhohats,[],1);
% vr22 = var(r2s,[],1);
% vgamma2 = var(gammas,[],1);
% 
% errorbar(lambdas,mrho2,sqrt(vrho2))
% % 
% x = lambdas'; y2 = mrho2'; err2 = sqrt(vrho2)';
% yplus2 = y2+err2; yminus2 = y2-err2;
% D = table(x,y2,err2,yplus2,yminus2);
% fName = sprintf('rho_N%i_T%i_beta%0.2f_seed%i.csv',N,T,beta,seed);
% writetable(D,fName)
% 
% y3 = predRhos';
% D3 = table(x,y3);
% fName = sprintf('predrho_N%i_T%i_beta%0.2f_seed%i.csv',N,T,beta,seed);
% writetable(D3,fName)



% 
% 
% x = lambdas'; y3 = mgammas'; err3 = sqrt(vgamma)';
% yplus3 = y3+err3; yminus3 = y3-err3;
% D = table(x,y3,err3,yplus3,yminus3);
% writetable(D,'REMbyLambda3.csv')
% 
% x = lambdas'; y4 = mgammas2'; err4 = sqrt(vgamma2)';
% yplus4 = y4+err4; yminus4 = y4-err4;
% D = table(x,y4,err4,yplus4,yminus4);
% writetable(D,'REMbyLambda4.csv')

% figure
% heatmap(Ns,betas,r2s)
% figure
% heatmap(Ns,betas,rhohats)
% figure
% heatmap(Ns,betas,rhos)

% figure
% heatmap(Ns,lambdas,r2s)
% figure
% heatmap(Ns,lambdas,rhohats)
% figure
% heatmap(Ns,lambdas,rhos)

% scatter(r2s,rhohats)
% hold on
% plot(r2s,rhohats)
% hold off


