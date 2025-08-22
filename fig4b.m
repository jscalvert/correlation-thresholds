%%%%%%%%%%%%%%%%%%%%%%%%%
%% PARTICLES ON A RING %%
%%%%%%%%%%%%%%%%%%%%%%%%%

figure

lambda = 0.5;


%% REFERENCE CURVE

L = 100;
x = 0:L;

qtilde = lambda.^min(x,L-x);

scatter(x./L,log(qtilde)./(L.*log(lambda)),40,'filled','o','MarkerFaceColor','k')

hold on

%% FIRST CURVE
L = 25;
x = 0:L;

p = 0.35;
color = p.*[48,182,243]./255 + (1-p).*[1,1,1];

% Set exit rates
q = lambda.^x + lambda.^(L-x);

% Fix the boundary rates
q(1) = lambda^L;
q(L+1) = lambda^L;

scatter(x./L,log(q)./(L.*log(lambda)),40,'filled','o','MarkerFaceColor',color,'MarkerFaceAlpha',1)

%% SECOND CURVE
L = 50;
x = 0:L;

p = 0.5;
color = p.*[48,182,243]./255 + (1-p).*[1,1,1];

% Set exit rates
q = lambda.^x + lambda.^(L-x);

% Fix the boundary rates
q(1) = lambda^L;
q(L+1) = lambda^L;

scatter(x./L,log(q)./(L.*log(lambda)),40,'filled','o','MarkerFaceColor',color,'MarkerFaceAlpha',1)

%% THIRD CURVE
L = 75;
x = 0:L

p = 1;
color = p.*[48,182,243]./255 + (1-p).*[1,1,1];

% Set exit rates
q = lambda.^x + lambda.^(L-x);

% Fix the boundary rates
q(1) = lambda^L;
q(L+1) = lambda^L;

scatter(x./L,log(q)./(L.*log(lambda)),40,'filled','o','MarkerFaceColor',color,'MarkerFaceAlpha',1)



set(gcf,'Color','white')
set(gca,'FontName','Times','FontSize',20,'LineWidth',1.5)
set(get(gca, 'XAxis'), 'FontWeight', 'bold');
set(get(gca, 'YAxis'), 'FontWeight', 'bold');
set(gca,'TickDir','in');
set(gca ,'Layer', 'Top');
set(gca,'TickLabelInterpreter','latex')
xlabel('$x/L$','FontSize',30,'Interpreter','latex')
ylabel('scaled exit rates','FontSize',30,'Interpreter','latex')
%ylim([0.5,1.1])
xticks([0,0.25,0.5,0.75,1])

set(gcf,'Position', [10 10 400 400])
hold off
%exportgraphics(gcf,'fig4b.pdf','ContentType','vector')
