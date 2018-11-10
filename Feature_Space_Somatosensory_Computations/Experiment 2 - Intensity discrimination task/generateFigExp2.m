clc; clear;
% Set some cool plot/graphics environments!
c =[ 0 0 1;...
     1 0 1;...
     1 0 0];
set(0,'defaultAxesColorOrder',c);
set(0,'defaultlinelinewidth',2)
set(0,'DefaultAxesFontSize',16)


% List of all subjects.
subject = {'1','2','3','4','5','6','7','8'};

% Generate .mat file extension for each subject.
sData = strcat('propcut_int_',subject,'.mat');

% Total number of subjects in the group analysis.
nSubject = length(subject);

% Use Gaussian CDF for fitting.
ft = fittype( '0.5*(1+erf((x-mu)/(sigma*sqrt(2))))', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( ft ); opts.Display = 'Off'; opts.Lower = [0 0]; opts.StartPoint = [0.4 0.02]; opts.Upper = [Inf Inf];

% Seven comparison intensities (ac) used in the experiment. 
% Standard intensity is always 0.4
ac = [0.05 0.20 0.32 0.40 0.48 0.60 0.75];

% Finer grid of the intensity range for smooth curve fitting done later.
acGrid = 0.05:0.005:0.75;

figure(1)

for k = 1:nSubject
    % Load kth subject's data.
    a = load(sData{k});
    for i = 1:7
        % Proabilities for each stimulus condition.
        % D1 and D2 only.
        p_d1(k,i) = mean(a.prob_d1{i}.prob);
        p_d2(k,i) = mean(a.prob_d2{i}.prob);
        % Near position.
        p_2(k,i) = mean(a.prob_near{1,i}.prob); 
        p_4(k,i) = mean(a.prob_near{2,i}.prob);
        p_7(k,i) = mean(a.prob_near{3,i}.prob);
        % Middle position.
        m_2(k,i) = mean(a.prob_middle{1,i}.prob); 
        m_4(k,i) = mean(a.prob_middle{2,i}.prob);
        m_7(k,i) = mean(a.prob_middle{3,i}.prob);
        % Far position.
        d_2(k,i) = mean(a.prob_far{1,i}.prob); 
        d_4(k,i) = mean(a.prob_far{2,i}.prob);
        d_7(k,i) = mean(a.prob_far{3,i}.prob);
    end
    
    % Plot fitted curve for each individual
    % Far
    subplot(231), hold on
    p_d1_ft = fit( [ac(1:3) ac(5:7)]', [p_d1(k,1:3) p_d1(k,5:7)]', ft, opts );
    plot(acGrid,p_d1_ft(acGrid),'k','LineWidth',1)
    d_2_ft = fit( ac', d_2(k,:)', ft, opts );
    plot(acGrid,d_2_ft(acGrid),'Color','b','LineWidth',1)
    d_4_ft = fit( ac', d_4(k,:)', ft, opts );
    plot(acGrid,d_4_ft(acGrid),'Color','m','LineWidth',1)
    d_7_ft = fit( ac', d_7(k,:)', ft, opts );
    plot(acGrid,d_7_ft(acGrid),'Color','r','LineWidth',1)
    
    % Middle
    subplot(232), hold on
    plot(acGrid,p_d1_ft(acGrid),'k','LineWidth',1)
    m_2_ft = fit( ac', m_2(k,:)', ft, opts );
    plot(acGrid,m_2_ft(acGrid),'Color','b','LineWidth',1)
    m_4_ft = fit( ac', m_4(k,:)', ft, opts );
    plot(acGrid,m_4_ft(acGrid),'Color','m','LineWidth',1)
    m_7_ft = fit( ac', m_7(k,:)', ft, opts );
    plot(acGrid,m_7_ft(acGrid),'Color','r','LineWidth',1)
    
    % Near
    subplot(233), hold on
    plot(acGrid,p_d1_ft(acGrid),'k','LineWidth',1)
    p_2_ft = fit( ac', p_2(k,:)', ft, opts );
    plot(acGrid,p_2_ft(acGrid),'Color','b','LineWidth',1)
    p_4_ft = fit( ac', p_4(k,:)', ft, opts );
    plot(acGrid,p_4_ft(acGrid),'Color','m','LineWidth',1)
    p_7_ft = fit( ac', p_7(k,:)', ft, opts );
    plot(acGrid,p_7_ft(acGrid),'Color','r','LineWidth',1)
end

% Mean probabilities (across subjects) for each stimulus condition.
prob_d1 = mean(p_d1);
near_2 = mean(p_2); near_4 = mean(p_4); near_7 = mean(p_7);
middle_2 = mean(m_2); middle_4 = mean(m_4); middle_7 = mean(m_7);
far_2 = mean(d_2); far_4 = mean(d_4); far_7 = mean(d_7);

% Starndard error of mean for each stimulus condition.
er_d1 = std(p_d1)/sqrt(nSubject);
erp_2 = std(p_2)/sqrt(nSubject); erp_4 = std(p_4)/sqrt(nSubject); erp_7 = std(p_7)/sqrt(nSubject);
erm_2 = std(m_2)/sqrt(nSubject); erm_4 = std(m_4)/sqrt(nSubject); erm_7 = std(m_7)/sqrt(nSubject);
erd_2 = std(d_2)/sqrt(nSubject); erd_4 = std(d_4)/sqrt(nSubject); erd_7 = std(d_7)/sqrt(nSubject);


% Fit mean probablities using the Gaussian CDF defined above.
% For D1 only, remove the data point when standard and comparison have same
% intensity (0.4). 
d1_fd_no = fit( [ac(1:3) ac(5:7)]', [prob_d1(1:3) prob_d1(5:7)]', ft, opts );

% Fit for near position.
near_fd_2 = fit( ac', near_2', ft, opts );
near_fd_4 = fit( ac', near_4', ft, opts );
near_fd_7 = fit( ac', near_7', ft, opts );

% Fit for middle position.
middle_fd_2 = fit( ac', middle_2', ft, opts );
middle_fd_4 = fit( ac', middle_4', ft, opts );
middle_fd_7 = fit( ac', middle_7', ft, opts );

% Fit for far position.
far_fd_2 = fit( ac', far_2', ft, opts );
far_fd_4 = fit( ac', far_4', ft, opts );
far_fd_7 = fit( ac', far_7', ft, opts );

% Plot fitted psychometric functions grouped by hand positions.
% Far position.
subplot(231), hold on
plot(acGrid,d1_fd_no(acGrid),'k')
h0 = errorbar([ac(1:3) ac(5:7)],[prob_d1(1:3) prob_d1(5:7)],[er_d1(1:3) er_d1(5:7)],'o'); set(h0,'Color','k','MarkerFaceColor','k','MarkerSize',1e-10,'LineWidth',2)
plot(acGrid,far_fd_2(acGrid),'Color','b')
h1 = errorbar(ac,far_2,erd_2,'o'); set(h1,'Color','b','MarkerFaceColor','b','MarkerSize',1e-10,'LineWidth',2)
plot(acGrid,far_fd_4(acGrid),'Color','m')
h2 = errorbar(ac,far_4,erd_4,'o'); set(h2,'Color','m','MarkerFaceColor','m','MarkerSize',1e-10,'LineWidth',2)
plot(acGrid,far_fd_7(acGrid),'Color','r')
h3 = errorbar(ac,far_7,erd_7,'o'); set(h3,'Color','r','MarkerFaceColor','r','MarkerSize',1e-10,'LineWidth',2)
errorbar_tick(h0,0); errorbar_tick(h1,0); errorbar_tick(h2,0); errorbar_tick(h3,0); 
axis([0.05 0.75 0 1]), set(gca,'YTick',[0.05 0.4 1],'XTick',[0.05 0.4 0.75])

xlabel('Comparison intensity (a.u.)')
ylabel({'P(Comparison', 'judged hihger)'})


% Middle position.
subplot(232), hold on
plot(acGrid,d1_fd_no(acGrid),'k')
h0 = errorbar([ac(1:3) ac(5:7)],[prob_d1(1:3) prob_d1(5:7)],[er_d1(1:3) er_d1(5:7)],'o'); set(h0,'Color','k','MarkerFaceColor','k','MarkerSize',1e-10,'LineWidth',2)
plot(acGrid,middle_fd_2(acGrid),'Color','b')
h1 = errorbar(ac,middle_2,erm_2,'o'); set(h1,'Color','b','MarkerFaceColor','b','MarkerSize',1e-10,'LineWidth',2)
plot(acGrid,middle_fd_4(acGrid),'Color','m')
h2 = errorbar(ac,middle_4,erm_4,'o'); set(h2,'Color','m','MarkerFaceColor','m','MarkerSize',1e-10,'LineWidth',2)
plot(acGrid,middle_fd_7(acGrid),'Color','r')
h3 = errorbar(ac,middle_7,erm_7,'o'); set(h3,'Color','r','MarkerFaceColor','r','MarkerSize',1e-10,'LineWidth',2)
errorbar_tick(h0,0); errorbar_tick(h1,0); errorbar_tick(h2,0); errorbar_tick(h3,0); 
axis([0.05 0.75 0 1]), set(gca,'YTick',[0.05 0.4 1],'XTick',[0.05 0.4 0.75])

% Near position.
subplot(233), hold on
plot(acGrid,d1_fd_no(acGrid),'k')
h0 = errorbar([ac(1:3) ac(5:7)],[prob_d1(1:3) prob_d1(5:7)],[er_d1(1:3) er_d1(5:7)],'o'); set(h0,'Color','k','MarkerFaceColor','k','MarkerSize',1e-10,'LineWidth',2)
plot(acGrid,near_fd_2(acGrid),'Color','b')
h1 = errorbar(ac,near_2,erp_2,'o'); set(h1,'Color','b','MarkerFaceColor','b','MarkerSize',1e-10,'LineWidth',2)
plot(acGrid,near_fd_4(acGrid),'Color','m')
h2 = errorbar(ac,near_4,erp_4,'o'); set(h2,'Color','m','MarkerFaceColor','m','MarkerSize',1e-10,'LineWidth',2)
plot(acGrid,near_fd_7(acGrid),'Color','r')
h3 = errorbar(ac,near_7,erp_7,'o'); set(h3,'Color','r','MarkerFaceColor','r','MarkerSize',1e-10,'LineWidth',2)
errorbar_tick(h0,0); errorbar_tick(h1,0); errorbar_tick(h2,0); errorbar_tick(h3,0); 
axis([0.05 0.75 0 1]), set(gca,'YTick',[0.05 0.4 1],'XTick',[0.05 0.4 0.75])

%%
% PSE and JND.
% Fit Gaussian CDF for all stimulus conditions for every subject.
% Extract PSE and JND from the fitted function.
for i = 1:nSubject
    % Fits for D1 and D2 only.
    D1only = fit( [ac(1:3) ac(5:7)]', p_d1(i,[1:3 5:7])', ft, opts );
    D2only = fit( [ac(1:3) ac(5:7)]', p_d2(i,[1:3 5:7])', ft, opts );
    % Fits for near position
    near_fd_2 = fit( ac', p_2(i,:)', ft, opts );
    near_fd_4 = fit( ac', p_4(i,:)', ft, opts );
    near_fd_7 = fit( ac', p_7(i,:)', ft, opts );
    % Fits for middle position.
    middle_fd_2 = fit( ac', m_2(i,:)', ft, opts );
    middle_fd_4 = fit( ac', m_4(i,:)', ft, opts );
    middle_fd_7 = fit( ac', m_7(i,:)', ft, opts );
    % Fits for far position.
    far_fd_2 = fit( ac', d_2(i,:)', ft, opts );
    far_fd_4 = fit( ac', d_4(i,:)', ft, opts );
    far_fd_7 = fit( ac', d_7(i,:)', ft, opts );
    % Extract PSE, indeed the change in PSE with respect to baseline PSE.
    mu_2(i,:) = ([far_fd_2.mu middle_fd_2.mu near_fd_2.mu] - D1only.mu);
    mu_4(i,:) = ([far_fd_4.mu middle_fd_4.mu near_fd_4.mu] - D1only.mu);
    mu_7(i,:) = ([far_fd_7.mu middle_fd_7.mu near_fd_7.mu] - D1only.mu);
    % Extract change in JND with respect to baseline JND.
    sigma_2(i,:) = ([far_fd_2.sigma middle_fd_2.sigma near_fd_2.sigma]-D1only.sigma);
    sigma_4(i,:) = ([far_fd_4.sigma middle_fd_4.sigma near_fd_4.sigma]-D1only.sigma);
    sigma_7(i,:) = ([far_fd_7.sigma middle_fd_7.sigma near_fd_7.sigma]-D1only.sigma);
end

% Combine PSE for each position; also the standard error of mean.
delMU = [mean(mu_2); mean(mu_4); mean(mu_7)];
semMU = [std(mu_2); std(mu_4); std(mu_7)]/sqrt(nSubject);
% Combine JND for each position; also the standard error of mean.
delSIGMA = [mean(sigma_2); mean(sigma_4); mean(sigma_7)];
semSIGMA = [std(sigma_2); std(sigma_4); std(sigma_7)]/sqrt(nSubject);


%% Plot PSE/Bias
clf;
% Group by hand locations
subplot(223), X = smartbar((delMU),0.9,1,{'Far','Middle','Near'},c); hold on
% Individual PSe
for subj = 1:8
    plot(X(1,:), (mu_2(subj,:)),'o','MarkerFaceColor','none','MarkerEdgeColor',[0.6 0.6 0.6])
    plot(X(2,:), (mu_4(subj,:)),'o','MarkerFaceColor','none','MarkerEdgeColor',[0.6 0.6 0.6])
    plot(X(3,:), (mu_7(subj,:)),'o','MarkerFaceColor','none','MarkerEdgeColor',[0.6 0.6 0.6])
end

% Error bars
eMU = errorbar(X(:),(delMU(:)),semMU(:),'.'); box off
errorbar_tick(eMU,25)
set(eMU,'Color','k','MarkerFaceColor','k','MarkerSize',1,'LineWidth',1.5)
ylabel('\DeltaPSE (Hz)')
ylim([-0.15 0.06])
set(gca,'YTick',[-0.15 0 0.06])

text(0,75,'100Hz','FontSize',14,'Color','b', 'HorizontalAlignment','left')
text(0,65,'200Hz','FontSize',14,'Color','m', 'HorizontalAlignment','left')
text(0,55,'300Hz','FontSize',14,'Color','r', 'HorizontalAlignment','left')


%% Plot of JND/Threshold
% Group by hand locations
subplot(224), X = smartbar(delSIGMA,0.9,1,{'Far','Middle','Near'},c); hold on
% Individual JND
for subj = 1:8
    plot(X(1,:), (sigma_2(subj,:)),'o','MarkerFaceColor','none','MarkerEdgeColor',[0.6 0.6 0.6])
    plot(X(2,:), (sigma_4(subj,:)),'o','MarkerFaceColor','none','MarkerEdgeColor',[0.6 0.6 0.6])
    plot(X(3,:), (sigma_7(subj,:)),'o','MarkerFaceColor','none','MarkerEdgeColor',[0.6 0.6 0.6])
end

% Error bars
eSIGMA = errorbar(X(:),delSIGMA(:),semSIGMA(:),'.'); box off
errorbar_tick(eSIGMA,25)
set(eSIGMA,'Color','k','MarkerFaceColor','k','MarkerSize',1,'LineWidth',1.5)
ylabel('\DeltaJND (Hz)')
ylim([-0.1 0.2])
set(gca,'YTick',[-0.1 0 0.1 0.2])
