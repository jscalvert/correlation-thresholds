%%%%%%%%%%%%%%%%%%%%%%%%%
%% PARTICLES ON A RING %%
%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Note that the use of corr in this script requires R2024a or later.
%
% Author: Jacob Calvert (calvert@gatech.edu)
% Date: February 10, 2025

% Set parameters
%Ls = [25,49,75,149]; % an odd integer >= 1
Ls = [50,75,100,125,150];
%Ls = [400];
%Ls = [25:5:150];
%Ls = 100;
numLs = numel(Ls);
% lambdas = exp([-4:0.00025:-0.00025,0.00025:0.00025:4]);
lambdas = exp([-1:0.0005:-0.0005,0.0005:0.0005:1]);
%lambdas = 1/1.1;
%lambdas = exp([-4:1:-1,1:1:4]);
%lambdas = exp([0.005:0.005:4]);
%lambdas = exp([-4:0.00125:-0.00125]);
%lambdas = 2;
%lambdas = [1.1];
numLambdas = numel(lambdas);

%color = [154,0,57];
%color = [255, 84, 0];
color = [57, 0, 153];
%color = [4,103,77];

rhos = zeros(numLs,numLambdas);
rhohats = zeros(numLs,numLambdas);
rs = zeros(numLs,numLambdas);

for i = 1:numLs
    L = Ls(i);
    % The states are 0 through L, which we denote by 1 through L+1
    x = 1:L+1;
    for j = 1:numLambdas
        lambda = lambdas(j);
    
        % Set exit rates
        q = zeros(1,L+1);
        p1 = zeros(1,L+1);
        p2 = zeros(1,L+1);
        for x = 1:L+1
            q(x) = lambda^(x-1) + lambda^(L-(x-1));
            p1(x) = (lambda^(x-1))/q(x); % left jump
            p2(x) = (lambda^(L-(x-1)))/q(x); % right jump
        end
        
        % Fix the boundaries
        q(1) = lambda^L;
        p1(1) = 0;
        p1(2) = 1;
        q(L+1) = lambda^L;
        p1(L+1) = 1;
        p2(L+1) = 0;
        
        % Determine the log stationary distribution
        %logp = zeros(1,L+1);
        
        z = 1:L+1;
        logp = log(lambda).*((z-1).*(L-(z-1)));
        logphatdiff = log(p2(1:end-1))-log(p1(2:end));
        logphat = [0,cumsum(logphatdiff)];
        
        % If p itself is needed, then:
        p = exp(logp)./sum(exp(logp));
        phat = exp(logphat)./sum(exp(logphat));
        
    
        % p = p./sum(p);
        % 
        % % The jump chain's stationary distribution is a function of p and q
        %phat = p.*q./sum(p.*q);
        % 
        % % Calculate variances
        % varlogphat = var(log(phat));
        % varlogq = var(log(q));
        % r = sqrt(varlogphat/varlogq)
        
        rhos(i,j) = corr(-logp',log(q)');
        rhohats(i,j) = corr(-logphat',log(q)');
        rs(i,j) = sqrt(var(logphat))/sqrt(var(log(q)));
        %rho2 = corr(-log(p)',log(q)')
        %rhohat = corrcoef(log(phat),-log(q))
        %i
    end
end
% 
% figure
% hold on
% for k = 1:numLs
%     colork = ((k/numLs).*color + (1 - k/numLs).*[255,255,255])./255;
%     plot(log(lambdas),rhos(k,:),'LineWidth',4,'Color',colork)
% end
% 
% plot(log(lambdas),(sqrt(15)/4).*ones(1,numLambdas),'LineWidth',2,'Color','k','LineStyle','--')
% plot(log(lambdas),(-sqrt(15)/4).*ones(1,numLambdas),'LineWidth',2,'Color','k','LineStyle','--')
% 
% set(gcf,'Color','white')
% set(gca,'FontName','Times','FontSize',20,'LineWidth',1.5)
% set(get(gca, 'XAxis'), 'FontWeight', 'bold');
% set(get(gca, 'YAxis'), 'FontWeight', 'bold');
% set(gca,'TickDir','in');
% set(gca ,'Layer', 'Top');
% set(gca,'TickLabelInterpreter','latex')
% xlabel('$\log \alpha$','FontSize',30,'Interpreter','latex')
% ylabel('$\rho$','FontSize',30,'Interpreter','latex')
% %xlim([-7,1])
% 
% set(gcf,'Position', [10 10 400 400])
% hold off

%xlim([])
%exportgraphics(gcf,'fig2b.pdf','ContentType','vector')




%figure
%scatter(log(q)',-logp')

%% NEW BLOCK
% figure
scatter(log(q)',-log(p)',50,'LineWidth',1.25,'MarkerEdgeColor',[51,84,122]./255,'Marker','o')
%scatter(log(q)',-log(p)',50,'LineWidth',1.25,'MarkerEdgeColor',color./255,'Marker','o')


set(gcf,'Color','white')
set(gca,'FontName','Times','FontSize',20,'LineWidth',1.5)
set(get(gca, 'XAxis'), 'FontWeight', 'bold');
set(get(gca, 'YAxis'), 'FontWeight', 'bold');
set(gca,'TickDir','in');
set(gca ,'Layer', 'Top');
set(gca,'TickLabelInterpreter','latex')
%xlabel('$\log q(x)$','FontSize',30,'Interpreter','latex')
%ylabel('$-\log \pi (x)$','FontSize',30,'Interpreter','latex')
xlabel('$\log q(x)$','FontSize',30,'Interpreter','latex')
ylabel('$-\log \pi (x)$','FontSize',30,'Interpreter','latex')
%xlim([-7,1])
xticks([-10:2:0])

set(gcf,'Position', [10 10 400 400])
%hold off
%hold on
% 
% exportgraphics(gcf,'fig2a.pdf','ContentType','vector')




% Calculate beta
%beta = 1 + rhohat*r;

%% Some plots

%figure
%plot(Ns',rhos)

%figure
% plot(Ls,rhohats,'black','LineWidth',2)
% hold on
% %subsamp = 1:1:numLambdas;
% %scatter(rs(subsamp),rhohats(subsamp),50,'black','LineWidth',1.25)
% scatter(Ls,rhohats,50,'red','LineWidth',1.25,'Marker','v')
% %ylim([-1,1])
% %xticks([0:2:10])
% %hold off

%exportgraphics(gcf,'fig2d.pdf','ContentType','vector')

%% Save results

% rhoData = [Ns',rhos];
% betaData = [Ns',betas];
% 
% if randomStateChoice == 1
%     csvwrite('FK_model_rho.csv',rhoData);
%     %csvwrite('FK_model_beta.csv',betaData);
% elseif randomStateChoice == 2
%     csvwrite('FK_model_rho_pi.csv',rhoData);
%     %csvwrite('FK_model_beta_pi.csv',betaData);
% end

