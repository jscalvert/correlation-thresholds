%% New wells and barriers

%% Set parameters
% seed
% 4 is nice
seed=4;
rng(seed)

% number of states
n = 10;
N = 2^n;

% magnitudes
sigmaW = 1;
sigmaB = 1;
sigmaF = 1;

%% Set adjacency matrix
% uncomment one of these

% % complete graph
%A = ones(N);
%A = A - diag(diag(A));

% % hypercube
A = double(hypercube(n));
A = A - diag(diag(A));

% % cycle on N states
%A = diag(ones(1,N-1),1)+diag(ones(1,N-1),-1);
%A(1,N) = 1; A(N,1) = 1;

%% Form energy wells and barriers
W = normrnd(0,sigmaW,1,N); % iid Gaussians with var = wellDepth
B = normrnd(0,sigmaB,N,N); % iid Gaussians with var = barrierHeight
F = normrnd(0,sigmaF,N,N); % iid Gaussians with var = forceStrength
B = triu(B,1) + (triu(B,1))'; % symmetric barriers
F = triu(F,1) - (triu(F,1))'; % antisymmetric forces
Q = zeros(N,N);
for i = 1:N
    Q(i,:) = exp(W(i) - B(i,:) + F(i,:)); % rates from wells, barriers, and forces
end

Q = Q.*A;

% Find exit rates
exitRates = sum(Q,2);
Q = Q - diag(exitRates);

% Find stationary distribution
statDist = pFromQ(Q);
statDist = statDist';

rho = corr(-log(statDist),log(exitRates))

scatter(log(exitRates),-log(statDist),75,'^','MarkerEdgeColor','#c59208','LineWidth',0.5,'MarkerEdgeAlpha',1)
hold on

% %% Make second plot
% 
% % magnitudes
% sigmaW = 1;
% sigmaB = 0;
% sigmaF = 4;
% 
% %% Form energy wells and barriers
% W = normrnd(0,sigmaW,1,N); % iid Gaussians with var = wellDepth
% B = normrnd(0,sigmaB,N,N); % iid Gaussians with var = barrierHeight
% F = normrnd(0,sigmaF,N,N); % iid Gaussians with var = forceStrength
% B = triu(B,1) + (triu(B,1))'; % symmetric barriers
% F = triu(F,1) - (triu(F,1))'; % antisymmetric forces
% Q = zeros(N,N);
% for i = 1:N
%     Q(i,:) = exp(W(i) - B(i,:) + F(i,:)); % rates from wells, barriers, and forces
% end
% 
% Q = Q.*A;
% 
% % Find exit rates
% exitRates = sum(Q,2);
% Q = Q - diag(exitRates);
% 
% % Find stationary distribution
% statDist = pFromQ(Q);
% statDist = statDist';
% 
% rho = corr(-log(statDist),log(exitRates))
% 
% scatter(log(exitRates),-log(statDist),75,'square','MarkerEdgeColor','#dd0152','LineWidth',0.5,'MarkerEdgeAlpha',1)

%% Make third plot
% Define transition rate matrix 
Q = zeros(N,N);
for i = 1:N
    Q(i,:) = exp(W(i));    
end

Q = Q.*A;

% Find exit rates
exitRates = sum(Q,2);
Q = Q - diag(exitRates);

% Find stationary distribution
% statDist = exp(-W)./sum(exp(-W));
statDist = pFromQ(Q);
statDist = statDist';

scatter(log(exitRates),-log(statDist),75,'o','MarkerEdgeColor','#000000','LineWidth',0.5,'MarkerEdgeAlpha',1)

set(gcf,'Color','white')
set(gca,'FontName','Times','FontSize',20,'LineWidth',1.5)
set(get(gca, 'XAxis'), 'FontWeight', 'bold');
set(get(gca, 'YAxis'), 'FontWeight', 'bold');
set(gca,'TickDir','in');
set(gca ,'Layer', 'Top');
set(gca,'TickLabelInterpreter','latex')
xlabel('$\log q_i$','FontSize',30,'Interpreter','latex')
ylabel('$-\log \pi_i$','FontSize',30,'Interpreter','latex')
%xlim([-7,1])

set(gcf,'Position', [10 10 400 400])
hold off
%xlim([])
exportgraphics(gcf,'fig1b.pdf','ContentType','vector')


%figure
%scatter(-log(exitRates),log(statDist),50,'square','MarkerEdgeColor',color1,'LineWidth',0.5,'MarkerEdgeAlpha',0.75)
%scatter(-log(exitRates),log(statDist),50,'diamond','MarkerEdgeColor',color2,'LineWidth',0.5,'MarkerEdgeAlpha',0.75)
%scatter(-log(exitRates),log(statDist),50,'o','MarkerEdgeColor',color3,'LineWidth',0.5,'MarkerEdgeAlpha',0.75)

%scatter(-log(exitRates),log(statDist),75,'square','MarkerEdgeColor','#000000','LineWidth',0.5,'MarkerEdgeAlpha',1)

%scatter(-log(exitRates),log(statDist),50,'o','MarkerEdgeColor','#3c0099','LineWidth',0.5,'MarkerEdgeAlpha',1)

%scatter(-log(exitRates),log(statDist),50,'diamond','MarkerEdgeColor','#8338ec','LineWidth',0.5,'MarkerEdgeAlpha',1)
%scatter(-log(exitRates),log(statDist),50,'o','MarkerEdgeColor','#3c0099','LineWidth',0.5,'MarkerEdgeAlpha',1)
%scatter(-log(exitRates),log(statDist),50,'diamond','MarkerEdgeColor','#c59208','LineWidth',0.5,'MarkerEdgeAlpha',1)

%scatter(-log(exitRates),log(statDist),50,'square','MarkerEdgeColor','#06946e','LineWidth',0.5,'MarkerEdgeAlpha',1)
