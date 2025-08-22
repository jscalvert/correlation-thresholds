%% New wells and barriers

%% Set parameters
% seed
% 4 is nice
seed=3;
rng(seed)

% Set key parameters
n = 10;
N = 2^n;
beta = 1;
h = 0;

T = 1;
lambdas = 0.5;
numLambdas = numel(lambdas);

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
end

scatter(log(exitRates),-log(statDist),50,'LineWidth',1,'MarkerEdgeColor',[51,84,122]./255,'Marker','o')

set(gcf,'Color','white')
set(gca,'FontName','Times','FontSize',20,'LineWidth',1.5)
set(get(gca, 'XAxis'), 'FontWeight', 'bold');
set(get(gca, 'YAxis'), 'FontWeight', 'bold');
set(gca,'TickDir','in');
set(gca ,'Layer', 'Top');
set(gca,'TickLabelInterpreter','latex')
xlabel('$\log q(x)$','FontSize',30,'Interpreter','latex')
ylabel('$-\log \pi (x)$','FontSize',30,'Interpreter','latex')
ylim([0,20])
set(gcf,'Position', [10 10 400 400])

%exportgraphics(gcf,'fig3a.pdf','ContentType','vector')