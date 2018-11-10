clc; clear;

% Set some cool plot/graphics environments!
c =[ 0 0 1;...
     1 0 0];
set(0,'defaultAxesColorOrder',c);
set(0,'defaultlinelinewidth',2)
set(0,'DefaultAxesFontSize',16)


subject = {'1','2','3','4','5','6','7','8'};

sData = strcat('propcut_freq_X',subject,'.mat');
n = 7; % number of comparison frequencies

unXed_d1 = zeros(1,n);
d1Xed = zeros(1,n);

unXed_100 = zeros(1,n); unXed_300 = zeros(1,n);
Xed_100 = zeros(1,n); Xed_300 = zeros(1,n);

nSubject = length(subject);
for k = 1:nSubject
    a = load(sData{k});
    for i = 1:7
        unXed_d1_all(k,i) = mean(a.prob_d1unXed{i}.prob);
        Xed_d1_all(k,i) = mean(a.prob_d1Xed{i}.prob);
        
        unXed_100_all(k,i) = mean(a.prob_unXed{1,i}.prob); % unXed
        unXed_300_all(k,i) = mean(a.prob_unXed{2,i}.prob);
        %
        Xed_100_all(k,i) = mean(a.prob_Xed{1,i}.prob); % Xed
        Xed_300_all(k,i) = mean(a.prob_Xed{2,i}.prob);
    end
end

unXed_d1 = mean(unXed_d1_all,1);
Xed_d1 = mean(Xed_d1_all,1);

unXed_100 = mean(unXed_100_all,1);
unXed_300 = mean(unXed_300_all,1);

Xed_100 = mean(Xed_100_all,1);
Xed_300 = mean(Xed_300_all,1);

sem_unXed_d1 = std(unXed_d1_all)/sqrt(nSubject);
sem_Xed_d1 = std(Xed_d1_all)/sqrt(nSubject);

sem_unXed_100 = std(unXed_100_all)/sqrt(nSubject); sem_unXed_300 = std(unXed_300_all)/sqrt(nSubject);
sem_Xed_100 = std(Xed_100_all)/sqrt(nSubject); sem_Xed_300 = std(Xed_300_all)/sqrt(nSubject);

ft = fittype( '0.5*(1+erf((x/1000-mu)/(sigma*sqrt(2))))', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( ft ); opts.Display = 'Off'; opts.Lower = [-Inf 0]; opts.StartPoint = [0.5 0.5]; opts.Upper = [Inf 1];

x = [100 140 180 200 220 260 300];
xp = 100:5:300;

unXed_ft_d1 = fit( [x(1:3) x(5:7)]', [unXed_d1(1:3) unXed_d1(5:7)]', ft, opts );
Xed_ft_d1 = fit( [x(1:3) x(5:7)]', [Xed_d1(1:3) Xed_d1(5:7)]', ft, opts );

unXed_ft_100 = fit( x', unXed_100', ft, opts );
unXed_ft_300 = fit( x', unXed_300', ft, opts );

Xed_ft_100 = fit( x', Xed_100', ft, opts );
Xed_ft_300 = fit( x', Xed_300', ft, opts );

clf; figure(1)


% Individual subject psychometric curve plots
for subj = 1:8
    unXed_myft_d1 = fit( [x(1:3) x(5:7)]', [unXed_d1_all(subj,1:3) unXed_d1_all(subj,5:7)]', ft, opts );
    Xed_myft_d1 = fit( [x(1:3) x(5:7)]', [Xed_d1_all(subj,1:3) Xed_d1_all(subj,5:7)]', ft, opts );
    
    unXed_myft_100 = fit( x', unXed_100_all(subj,:)', ft, opts );
    unXed_myft_300 = fit( x', unXed_300_all(subj,:)', ft, opts );
    
    Xed_myft_100 = fit( x', Xed_100_all(subj,:)', ft, opts );
    Xed_myft_300 = fit( x', Xed_300_all(subj,:)', ft, opts );
    
    subplot(221), hold on
    plot(xp,unXed_myft_d1(xp),'k','LineWidth',1)
    plot(xp,unXed_myft_100(xp),'Color',c(1,:),'LineWidth',1)
    plot(xp,unXed_myft_300(xp),'Color',c(2,:),'LineWidth',1)
    
    subplot(222), hold on
    plot(xp,Xed_myft_d1(xp),'k','LineWidth',1)
    plot(xp,Xed_myft_100(xp),'Color',c(1,:),'LineWidth',1)
    plot(xp,Xed_myft_300(xp),'Color',c(2,:),'LineWidth',1)
end





subplot(221)
% plot(x,unXed_100,'o', x,unXed_300,'o', x([1:3 5:7]),unXed_d1([1:3 5:7]),'ok'), hold on
plot(xp,unXed_ft_100(xp), xp,unXed_ft_300(xp)), hold on
plot(xp,unXed_ft_d1(xp),'k')

if nSubject > 1
    h1 = errorbar(x,unXed_100,sem_unXed_100,'o'); set(h1,'Color',c(1,:),'MarkerFaceColor',c(1,:),'MarkerSize',1e-10,'LineWidth',1.5)
    h3 = errorbar(x,unXed_300,sem_unXed_300,'o'); set(h3,'Color',c(2,:),'MarkerFaceColor',c(2,:),'MarkerSize',1e-10,'LineWidth',1.5)
    h0 = errorbar([x(1:3) x(5:7)],[unXed_d1(1:3) unXed_d1(5:7)],[sem_unXed_d1(1:3) sem_unXed_d1(5:7)],'o'); set(h0,'Color','k','MarkerFaceColor','k','MarkerSize',1e-10,'LineWidth',1.5)
    errorbar_tick(h1,0)
    errorbar_tick(h3,0)
    errorbar_tick(h0,0)
end
axis([100 300 -0.08 1]), set(gca,'YTick',[0 0.5 1])
xlabel('Comparison frequency (Hz)'), ylabel('P(judged Fc > Fs)')
title('Uncrossed'), box off

text(100,1.00,'Fd = 100 Hz','FontSize',15,'Color',c(1,:), 'HorizontalAlignment','left')
text(100,0.85,'Fd = 300 Hz','FontSize',15,'Color',c(2,:), 'HorizontalAlignment','left')
text(100,0.70,'Baseline','FontSize',15,'Color','k', 'HorizontalAlignment','left')

text(305,0.06,'Fs = 200 Hz','FontSize',15,'Color','k', 'HorizontalAlignment','right')


subplot(222)
% plot(x,Xed_100,'o', x,Xed_300,'o', x([1:3 5:7]),Xed_d1([1:3 5:7]),'ok'), hold on
plot(xp,Xed_ft_100(xp), xp,Xed_ft_300(xp)), hold on
plot(xp,Xed_ft_d1(xp),'k')

if nSubject > 1
    h1 = errorbar(x,Xed_100,sem_Xed_100,'o'); set(h1,'Color',c(1,:),'MarkerFaceColor',c(1,:),'MarkerSize',1e-10,'LineWidth',1.5)
    h3 = errorbar(x,Xed_300,sem_Xed_300,'o'); set(h3,'Color',c(2,:),'MarkerFaceColor',c(2,:),'MarkerSize',1e-10,'LineWidth',1.5)
    h0 = errorbar([x(1:3) x(5:7)],[Xed_d1(1:3) Xed_d1(5:7)],[sem_Xed_d1(1:3) sem_Xed_d1(5:7)],'o'); set(h0,'Color','k','MarkerFaceColor','k','MarkerSize',1e-10,'LineWidth',1.5)
    errorbar_tick(h1,0)
    errorbar_tick(h3,0)
    errorbar_tick(h0,0)
end
axis([100 300 -0.08 1]), set(gca,'YTick',[0 0.5 1])
xlabel('Comparison frequency (Hz)'), ylabel('P(judged Fc > Fs)')
title('Crossed'), box off

for i = 1:nSubject
    unXed_ft_d1 = fit( [x(1:3) x(5:7)]', unXed_d1_all(i,[1:3 5:7])', ft, opts );
    Xed_ft_d1 = fit( [x(1:3) x(5:7)]', Xed_d1_all(i,[1:3 5:7])', ft, opts );
    
    unXed_ft_100 = fit( x', unXed_100_all(i,:)', ft, opts );
    unXed_ft_300 = fit( x', unXed_300_all(i,:)', ft, opts );
    
    Xed_ft_100 = fit( x', Xed_100_all(i,:)', ft, opts );
    Xed_ft_300 = fit( x', Xed_300_all(i,:)', ft, opts );
    
    pse_unXed(i,:) = 1000 * ([unXed_ft_100.mu unXed_ft_300.mu unXed_ft_d1.mu]);
    pse_Xed(i,:) = 1000 * ([Xed_ft_100.mu Xed_ft_300.mu  Xed_ft_d1.mu]);
    
    jnd_unXed(i,:) = 1000 * ([unXed_ft_100.sigma unXed_ft_300.sigma unXed_ft_d1.sigma]);
    jnd_Xed(i,:) = 1000 * ([Xed_ft_100.sigma Xed_ft_300.sigma Xed_ft_d1.sigma]);
end

% [(1:nSubject)' pse_Xed pse_unXed]
% [(1:nSubject)' jnd_Xed jnd_unXed]

dpseXed = bsxfun(@minus,pse_Xed(:,1:2),pse_Xed(:,end));
dpseunXed = bsxfun(@minus,pse_unXed(:,1:2),pse_unXed(:,end));

djndXed = bsxfun(@minus,jnd_Xed(:,1:2),jnd_Xed(:,end));
djndunXed = bsxfun(@minus,jnd_unXed(:,1:2),jnd_unXed(:,end));

% [(1:nSubject)' dpseXed dpseunXed]
% [(1:nSubject)' djndXed djndunXed]

avg_del_pse_unXed = mean(dpseunXed, 1);
avg_del_pse_Xed = mean(dpseXed, 1);
delPSE = [avg_del_pse_unXed' avg_del_pse_Xed'];
semPSE = [std(dpseunXed)' std(dpseXed)']/sqrt(nSubject);

avg_del_jnd_unXed = mean(djndunXed,1);
avg_del_jnd_Xed = mean(djndXed,1);
delJND = [avg_del_jnd_unXed' avg_del_jnd_Xed'];
semJND = [std(djndunXed)' std(djndXed)']/sqrt(nSubject);


subplot(223), X = smartbar((delPSE),0.9,1,{'Uncrossed','Crossed'},c); hold on
if nSubject > 1
    ebPSE = errorbar(X(:),(delPSE(:)),semPSE(:),'o');
    set(ebPSE,'Color','k','MarkerFaceColor','k','MarkerSize',4,'LineWidth',1.5)
    errorbar_tick(ebPSE,20)
end

for subj = 1:8
    plot(X(:,1), [dpseunXed(subj,1) dpseunXed(subj,2)],'o','MarkerFaceColor','none','MarkerEdgeColor',[0.6 0.6 0.6])
    plot(X(:,2),[dpseXed(subj,1) dpseXed(subj,2)],'o','MarkerFaceColor','none','MarkerEdgeColor',[0.6 0.6 0.6])
end

% xlabel('Distractor (Hz)')
ylabel('\DeltaPSE (Hz)')
ylim([-60 60]); box off
set(gca,'YTick',[-60 0 60])
% text(0,40,'100 Hz','FontSize',14,'Color',c(1,:), 'HorizontalAlignment','left')
% text(1,40,'300 Hz','FontSize',14,'Color',c(2,:), 'HorizontalAlignment','left')
% text(4,50,'Standard = 200 Hz','FontSize',13, 'HorizontalAlignment','right')


subplot(224), X = smartbar(delJND,0.9,1,{'Uncrossed','Crossed'},c); hold on
if nSubject > 1
    ebJND = errorbar(X(:),delJND(:),semJND(:),'o');
    set(ebJND,'Color','k','MarkerFaceColor','k','MarkerSize',4,'LineWidth',1.5)
    errorbar_tick(ebJND,20)
end

for subj = 1:8
    plot(X(:,1), [djndunXed(subj,1) djndunXed(subj,2)],'o','MarkerFaceColor','none','MarkerEdgeColor',[0.6 0.6 0.6])
    plot(X(:,2),[djndXed(subj,1) djndXed(subj,2)],'o','MarkerFaceColor','none','MarkerEdgeColor',[0.6 0.6 0.6])
end


ylabel('\DeltaJND (Hz)')
axis([-0.6 4.5 -20 70]); box off
set(gca,'YTick',[0 70])

% text(4.5,40,'Distractor','FontSize',15,'Color','k', 'HorizontalAlignment','right')
% text(4.5,34,'100 Hz','FontSize',15,'Color',c(1,:), 'HorizontalAlignment','right')
% text(4.5,28,'300 Hz','FontSize',15,'Color',c(2,:), 'HorizontalAlignment','right')





