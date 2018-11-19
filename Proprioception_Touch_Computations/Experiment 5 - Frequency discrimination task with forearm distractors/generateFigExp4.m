clc; clear;

% Set some cool plot/graphics environments!
c =[ 0 0 1;...
     1 0 0];
set(0,'defaultAxesColorOrder',c);
set(0,'defaultlinelinewidth',2)
set(0,'DefaultAxesFontSize',16)

subject = {'1','2','3','4','5','6','7','8'};

sData = strcat('propcut_freq_digitForearm',subject,'.mat');
n = 7; % number of comparison frequencies

d1only = zeros(1,n);

far_100 = zeros(1,n); far_300 = zeros(1,n);
near_100 = zeros(1,n); near_300 = zeros(1,n);

nSubject = length(subject);
for k = 1:nSubject
    a = load(sData{k});
    for i = 1:7
        d1only_all(k,i) = mean(a.prob_d1only{i}.prob);
        
        far_100_all(k,i) = mean(a.prob_far{1,i}.prob); % far
        far_300_all(k,i) = mean(a.prob_far{2,i}.prob);
        %
        near_100_all(k,i) = mean(a.prob_near{1,i}.prob); % near
        near_300_all(k,i) = mean(a.prob_near{2,i}.prob);
    end
end

d1only = mean(d1only_all,1); 

far_100 = mean(far_100_all,1);
far_300 = mean(far_300_all,1);

near_100 = mean(near_100_all,1);
near_300 = mean(near_300_all,1);

sem_d1only = std(d1only_all)/sqrt(nSubject);

sem_far_100 = std(far_100_all)/sqrt(nSubject); sem_far_300 = std(far_300_all)/sqrt(nSubject);
sem_near_100 = std(near_100_all)/sqrt(nSubject); sem_near_300 = std(near_300_all)/sqrt(nSubject);

ft = fittype( '0.5*(1+erf((x/1000-mu)/(sigma*sqrt(2))))', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( ft ); opts.Display = 'Off'; opts.Lower = [-Inf 0]; opts.StartPoint = [0.5 0.5]; opts.Upper = [Inf 1];

x = [100 140 180 200 220 260 300];
xp = 100:10:300;

d1only_ft = fit( [x(1:3) x(5:7)]', [d1only(1:3) d1only(5:7)]', ft, opts );

far_ft_100 = fit( x', far_100', ft, opts );
far_ft_300 = fit( x', far_300', ft, opts );

near_ft_100 = fit( x', near_100', ft, opts );
near_ft_300 = fit( x', near_300', ft, opts );

clf; figure(1)

for subj = 1:8
    d1only_myft = fit( [x(1:3) x(5:7)]', [d1only_all(subj,1:3) d1only_all(subj,5:7)]', ft, opts );
    
    far_myft_100 = fit( x', far_100_all(subj,:)', ft, opts );
    far_myft_300 = fit( x', far_300_all(subj,:)', ft, opts );
    
    near_myft_100 = fit( x', near_100_all(subj,:)', ft, opts );
    near_myft_300 = fit( x', near_300_all(subj,:)', ft, opts );
    
    subplot(223)
    plot(xp, d1only_myft(xp), 'k', 'LineWidth', 1), hold on
    plot(xp, near_myft_100(xp), 'Color',c(1,:), 'LineWidth', 1)
    plot(xp, near_myft_300(xp), 'Color',c(2,:), 'LineWidth', 1)
    
    subplot(221)
    plot(xp, d1only_myft(xp), 'k', 'LineWidth', 1), hold on
    plot(xp, far_myft_100(xp), 'Color',c(1,:), 'LineWidth', 1)
    plot(xp, far_myft_300(xp), 'Color',c(2,:), 'LineWidth', 1)
    
end


subplot(2,2,3)
plot(xp,near_ft_100(xp), xp,near_ft_300(xp),...
    xp,d1only_ft(xp),'k'), hold on
% plot(x,near_100,'ob','MarkerFaceColor','b');
% plot(x,near_300,'or','MarkerFaceColor','r');
plot(x([1:3 5:7]),d1only([1:3 5:7]),'ok','MarkerFaceColor','k');


if nSubject > 1
    h100near = errorbar(x,near_100,sem_near_100,'o'); set(h100near,'Color',c(1,:),'MarkerFaceColor',c(1,:),'MarkerSize',1e-10,'LineWidth',1.5); errorbar_tick(h100near,70)
    h300near = errorbar(x,near_300,sem_near_300,'o'); set(h300near,'Color',c(2,:),'MarkerFaceColor',c(2,:),'MarkerSize',1e-10,'LineWidth',1.5); errorbar_tick(h300near,70)
    h0 = errorbar([x(1:3) x(5:7)],[d1only(1:3) d1only(5:7)],[sem_d1only(1:3) sem_d1only(5:7)],'o'); set(h0,'Color','k','MarkerFaceColor','k','MarkerSize',1e-10,'LineWidth',1.5)
    errorbar_tick(h100near,0)
    errorbar_tick(h300near,0)
    errorbar_tick(h0,0)
end
axis([100 300 -0.08 1]), set(gca,'YTick',[0 0.5 1])
xlabel('Fc (Hz)'), ylabel('P(judged Fc > Fs)'), box off
title('Near')



subplot(2,2,1)
plot(xp,far_ft_100(xp), xp,far_ft_300(xp),...
    xp,d1only_ft(xp),'k'), hold on
% plot(x,far_100,'ob', x,far_300,'or')

if nSubject > 1
    h100far = errorbar(x,far_100,sem_far_100,'o'); set(h100far,'Color',c(1,:),'MarkerFaceColor',c(1,:),'MarkerSize',1e-10,'LineWidth',1.5); errorbar_tick(h100far,70)
    h300far = errorbar(x,far_300,sem_far_300,'o'); set(h300far,'Color',c(2,:),'MarkerFaceColor',c(2,:),'MarkerSize',1e-10,'LineWidth',1.5); errorbar_tick(h300far,70)
    h0 = errorbar([x(1:3) x(5:7)],[d1only(1:3) d1only(5:7)],[sem_d1only(1:3) sem_d1only(5:7)],'o'); set(h0,'Color','k','MarkerFaceColor','k','MarkerSize',1e-10,'LineWidth',1.5)
    errorbar_tick(h100far,0)
    errorbar_tick(h300far,0)
    errorbar_tick(h0,0)
end
axis([100 300 -0.08 1]), set(gca,'YTick',[0 0.5 1])
xlabel('Fc (Hz)'), ylabel('P(judged Fc > Fs)'), box off
title('Far')

text(100,1.00,'Fd = 100 Hz','FontSize',15,'Color',c(1,:), 'HorizontalAlignment','left')
text(100,0.85,'Fd = 300 Hz','FontSize',15,'Color',c(2,:), 'HorizontalAlignment','left')
text(100,0.70,'Baseline','FontSize',15,'Color','k', 'HorizontalAlignment','left')
text(305,0.06,'Fs = 200 Hz','FontSize',15,'Color','k', 'HorizontalAlignment','right')



for i = 1:nSubject
    d1only_ft = fit( [x(1:3) x(5:7)]', d1only_all(i,[1:3 5:7])', ft, opts );
    
    far_ft_100 = fit( x', far_100_all(i,:)', ft, opts );
    far_ft_300 = fit( x', far_300_all(i,:)', ft, opts );
    
    near_ft_100 = fit( x', near_100_all(i,:)', ft, opts );
    near_ft_300 = fit( x', near_300_all(i,:)', ft, opts );
    
    pse_far(i,:) = 1000 * ([far_ft_100.mu far_ft_300.mu d1only_ft.mu]);
    pse_near(i,:) = 1000 * ([near_ft_100.mu near_ft_300.mu  d1only_ft.mu]);
    
    jnd_far(i,:) = 1000 * ([far_ft_100.sigma far_ft_300.sigma d1only_ft.sigma]);
    jnd_near(i,:) = 1000 * ([near_ft_100.sigma near_ft_300.sigma d1only_ft.sigma]);
end


dpsenear = bsxfun(@minus,pse_near(:,1:2),pse_near(:,end));
dpsefar = bsxfun(@minus,pse_far(:,1:2),pse_far(:,end));

djndnear = bsxfun(@minus,jnd_near(:,1:2),jnd_near(:,end));
djndfar = bsxfun(@minus,jnd_far(:,1:2),jnd_far(:,end));


avg_del_pse_far = mean(dpsefar, 1);
avg_del_pse_near = mean(dpsenear, 1);
delPSE = [avg_del_pse_far' avg_del_pse_near'];
semPSE = [std(dpsefar)' std(dpsenear)']/sqrt(nSubject);

avg_del_jnd_near = mean(djndnear,1);
avg_del_jnd_far = mean(djndfar,1);
delJND = [avg_del_jnd_far' avg_del_jnd_near'];
semJND = [std(djndfar)' std(djndnear)']/sqrt(nSubject);


subplot(222), X = smartbar((delPSE),0.9,1,{'Far','Near'},c); hold on
if nSubject > 1
    ebPSE = errorbar(X(:),(delPSE(:)),semPSE(:),'o');
    set(ebPSE,'Color','k','MarkerFaceColor','k','MarkerSize',4,'LineWidth',1.5)
    errorbar_tick(ebPSE,20)
end

for subj = 1:8
    plot(X(:,1), [dpsefar(subj,1) dpsefar(subj,2)],'o','MarkerFaceColor','none','MarkerEdgeColor',[0.6 0.6 0.6])
    plot(X(:,2),[dpsenear(subj,1) dpsenear(subj,2)],'o','MarkerFaceColor','none','MarkerEdgeColor',[0.6 0.6 0.6])
end

% xlabel('Distractor (Hz)'), 
ylabel('\DeltaPSE (Hz)')
ylim([-60 60]); box off
set(gca,'YTick',[-60 0 60])
% text(0,40,'100 Hz','FontSize',14,'Color',c(1,:), 'HorizontalAlignment','left')
% text(1,40,'300 Hz','FontSize',14,'Color',c(2,:), 'HorizontalAlignment','left')
% text(4,50,'Standard = 200 Hz','FontSize',13, 'HorizontalAlignment','right')


subplot(224), X = smartbar(delJND,0.9,1,{'Far','Near'},c); hold on
if nSubject > 1
    ebJND = errorbar(X(:),delJND(:),semJND(:),'o');
    set(ebJND,'Color','k','MarkerFaceColor','k','MarkerSize',4,'LineWidth',1.5)
    errorbar_tick(ebJND,20)
end

for subj = 1:8
    plot(X(:,1), [djndfar(subj,1) djndfar(subj,2)],'o','MarkerFaceColor','none','MarkerEdgeColor',[0.6 0.6 0.6])
    plot(X(:,2),[djndnear(subj,1) djndnear(subj,2)],'o','MarkerFaceColor','none','MarkerEdgeColor',[0.6 0.6 0.6])
end

% xlabel('Distractor'), 
ylabel('\DeltaJND (Hz)')
axis([-0.6 4.5 0 50]); box off
set(gca,'YTick',[0 50])

