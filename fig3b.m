%% New wells and barriers

%% Set parameters
% seed
% 4 is nice
seed=3;
rng(seed)

% Set key parameters
n = 4;
N = 2^n;
beta = 1;
h = 0;

T = 50;
lambdas = [0:0.1:1];
numLambdas = numel(lambdas);

color = [167,214,78]./255;

rhos = zeros(T,numel(lambdas));
rhohats = zeros(T,numel(lambdas));
rs = zeros(T,numel(lambdas));

%% Set adjacency matrix

% % hypercube
A = double(hypercube(n));
A = A - diag(diag(A));

for k = 1:numLambdas
    lambda = lambdas(k);
    for t = 1:T
        %% Form energy wells and barriers
        % Couplings
        J = normrnd(0,1/sqrt(n),n,n);
        
        % Energies
        E = SKEnergy(n,J,h);
        % Rates
        Q = zeros(N);
        for l=1:N
            for m=1:N
                Q(l,m) = exp(beta*(lambda*(E(l)+E(m))-E(m)));
            end
        end
        
        Q = Q.*A;
        
        % Find exit rates
        exitRates = sum(Q,2);
        Q = Q - diag(exitRates);
        
        % Find stationary distribution
        statDist = pFromQ(Q);
        statDist = statDist';
        phat = (statDist.*exitRates)./sum(statDist.*exitRates);

        rhos(t,k) = corr(-log(statDist),log(exitRates));
        rhohats(t,k) = corr(-log(phat),log(exitRates));
        rs(t,k) = sqrt(var(log(phat))/var(log(exitRates)));
    end
    k
end

figure
hold on
% for k = 1:numLs
%     colork = ((k/numLs).*color + (1 - k/numLs).*[255,255,255])./255;
%     plot(log(lambdas),rhos(k,:),'LineWidth',4,'Color',colork)
% end
mrho = mean(rhos);
vrho = var(rhos);
errorbar(lambdas,mrho,sqrt(vrho),sqrt(vrho),'LineWidth',4,'Color',color)

%plot(lambdas,mrho,'LineWidth',4,'Color',color)


predRhos = sign(2.*lambdas-1 + 4.*(1-lambdas)./n);
plot(lambdas,predRhos,'LineWidth',3,'Color',color,'LineStyle','--')

%% MAKE SECOND PLOT

% Set key parameters
n = 10;
N = 2^n;
beta = 1;
h = 0;
T = 10;

lambdas = [0:0.1:1];
numLambdas = numel(lambdas);

color = [57, 0, 153]./255;

rhos = zeros(T,numel(lambdas));
rhohats = zeros(T,numel(lambdas));
rs = zeros(T,numel(lambdas));

%% Set adjacency matrix

% % hypercube
A = double(hypercube(n));
A = A - diag(diag(A));

for k = 1:numLambdas
    lambda = lambdas(k);
    for t = 1:T
        %% Form energy wells and barriers
        % Couplings
        J = normrnd(0,1/sqrt(n),n,n);

        % Energies
        E = SKEnergy(n,J,h);
        % Rates
        Q = zeros(N);
        for l=1:N
            for m=1:N
                Q(l,m) = exp(beta*(lambda*(E(l)+E(m))-E(m)));
            end
        end

        Q = Q.*A;

        % Find exit rates
        exitRates = sum(Q,2);
        Q = Q - diag(exitRates);

        % Find stationary distribution
        statDist = pFromQ(Q);
        statDist = statDist';
        phat = (statDist.*exitRates)./sum(statDist.*exitRates);

        rhos(t,k) = corr(-log(statDist),log(exitRates));
        rhohats(t,k) = corr(-log(phat),log(exitRates));
        rs(t,k) = sqrt(var(log(phat))/var(log(exitRates)));
    end
    k
end

%figure
%hold on
% for k = 1:numLs
%     colork = ((k/numLs).*color + (1 - k/numLs).*[255,255,255])./255;
%     plot(log(lambdas),rhos(k,:),'LineWidth',4,'Color',colork)
% end
mrho = mean(rhos);
vrho = var(rhos);
errorbar(lambdas,mrho,sqrt(vrho),sqrt(vrho),'LineWidth',4,'Color',color)

%plot(lambdas,mrho,'LineWidth',4,'Color',color)

predRhos = sign(2.*lambdas-1 + 4.*(1-lambdas)./n);
plot(lambdas,predRhos,'LineWidth',3,'Color',color,'LineStyle','--')

%% UPDATE PLOT


set(gcf,'Color','white')
set(gca,'FontName','Times','FontSize',20,'LineWidth',1.5)
set(get(gca, 'XAxis'), 'FontWeight', 'bold');
set(get(gca, 'YAxis'), 'FontWeight', 'bold');
set(gca,'TickDir','in');
set(gca ,'Layer', 'Top');
set(gca,'TickLabelInterpreter','latex')
xlabel('$\lambda$','FontSize',30,'Interpreter','latex')
ylabel('$\rho$','FontSize',30,'Interpreter','latex','Rotation',0)
%xticks([0,0.25,0.5,0.75,1])

%ylim([0,20])
set(gcf,'Position', [10 10 400 400])
ylim([-1,1])

exportgraphics(gcf,'fig3b_v3.pdf','ContentType','vector')