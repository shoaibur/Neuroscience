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
sData = strcat('propcut_freq_',subject,'.mat');

% Total number of subjects in the group analysis.
nSubject = length(subject);

% Use Gaussian CDF for fitting.
ft = fittype( '0.5*(1+erf((x/1000-mu)/(sigma*sqrt(2))))', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( ft ); opts.Display = 'Off'; opts.Lower = [-Inf 0]; opts.StartPoint = [0.5 0.5]; opts.Upper = [Inf 1];

% Seven comparison frequenicies (fc) used in the experiment. 
% Standard frequency is always 200 Hz.
fc = [100 140 180 200 220 260 300];

% Finer grid of the fequency range for smooth curve fitting done later.
fcGrid = 100:5:300;

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
        p_100(k,i) = mean(a.prob_near{1,i}.prob); 
        p_200(k,i) = mean(a.prob_near{2,i}.prob);
        p_300(k,i) = mean(a.prob_near{3,i}.prob);
        % Middle position.
        m_100(k,i) = mean(a.prob_middle{1,i}.prob); 
        m_200(k,i) = mean(a.prob_middle{2,i}.prob);
        m_300(k,i) = mean(a.prob_middle{3,i}.prob);
        % Far position.
        d_100(k,i) = mean(a.prob_far{1,i}.prob); 
        d_200(k,i) = mean(a.prob_far{2,i}.prob);
        d_300(k,i) = mean(a.prob_far{3,i}.prob);
    end
    
    % Plot fitted curve for each individual
    % Far
    subplot(231), hold on
    p_d1_ft = fit( [fc(1:3) fc(5:7)]', [p_d1(k,1:3) p_d1(k,5:7)]', ft, opts );
    plot(fcGrid,p_d1_ft(fcGrid),'k','LineWidth',1)
    d_100_ft = fit( fc', d_100(k,:)', ft, opts );
    plot(fcGrid,d_100_ft(fcGrid),'Color','b','LineWidth',1)
    d_200_ft = fit( fc', d_200(k,:)', ft, opts );
    plot(fcGrid,d_200_ft(fcGrid),'Color','m','LineWidth',1)
    d_300_ft = fit( fc', d_300(k,:)', ft, opts );
    plot(fcGrid,d_300_ft(fcGrid),'Color','r','LineWidth',1)
    
    % Middle
    subplot(232), hold on
    plot(fcGrid,p_d1_ft(fcGrid),'k','LineWidth',1)
    m_100_ft = fit( fc', m_100(k,:)', ft, opts );
    plot(fcGrid,m_100_ft(fcGrid),'Color','b','LineWidth',1)
    m_200_ft = fit( fc', m_200(k,:)', ft, opts );
    plot(fcGrid,m_200_ft(fcGrid),'Color','m','LineWidth',1)
    m_300_ft = fit( fc', m_300(k,:)', ft, opts );
    plot(fcGrid,m_300_ft(fcGrid),'Color','r','LineWidth',1)
    
    % Near
    subplot(233), hold on
    plot(fcGrid,p_d1_ft(fcGrid),'k','LineWidth',1)
    p_100_ft = fit( fc', p_100(k,:)', ft, opts );
    plot(fcGrid,p_100_ft(fcGrid),'Color','b','LineWidth',1)
    p_200_ft = fit( fc', p_200(k,:)', ft, opts );
    plot(fcGrid,p_200_ft(fcGrid),'Color','m','LineWidth',1)
    p_300_ft = fit( fc', p_300(k,:)', ft, opts );
    plot(fcGrid,p_300_ft(fcGrid),'Color','r','LineWidth',1)
end

% Mean probabilities (across subjects) for each stimulus condition.
prob_d1 = mean(p_d1);
near_100 = mean(p_100); near_200 = mean(p_200); near_300 = mean(p_300);
middle_100 = mean(m_100); middle_200 = mean(m_200); middle_300 = mean(m_300);
far_100 = mean(d_100); far_200 = mean(d_200); far_300 = mean(d_300);

% Starndard error of mean for each stimulus condition.
er_d1 = std(p_d1)/sqrt(nSubject);
erp_100 = std(p_100)/sqrt(nSubject); erp_200 = std(p_200)/sqrt(nSubject); erp_300 = std(p_300)/sqrt(nSubject);
erm_100 = std(m_100)/sqrt(nSubject); erm_200 = std(m_200)/sqrt(nSubject); erm_300 = std(m_300)/sqrt(nSubject);
erd_100 = std(d_100)/sqrt(nSubject); erd_200 = std(d_200)/sqrt(nSubject); erd_300 = std(d_300)/sqrt(nSubject);


% Fit mean probablities using the Gaussian CDF defined above.
% For D1 only, remove the data point when standard and comparison have same
% frequency (200 Hz). 
d1_fd_no = fit( [fc(1:3) fc(5:7)]', [prob_d1(1:3) prob_d1(5:7)]', ft, opts );

% Fit for near position.
near_fd_100 = fit( fc', near_100', ft, opts );
near_fd_200 = fit( fc', near_200', ft, opts );
near_fd_300 = fit( fc', near_300', ft, opts );

% Fit for middle position.
middle_fd_100 = fit( fc', middle_100', ft, opts );
middle_fd_200 = fit( fc', middle_200', ft, opts );
middle_fd_300 = fit( fc', middle_300', ft, opts );

% Fit for far position.
far_fd_100 = fit( fc', far_100', ft, opts );
far_fd_200 = fit( fc', far_200', ft, opts );
far_fd_300 = fit( fc', far_300', ft, opts );

% Plot fitted psychometric functions grouped by hand positions.
% Far position.
subplot(231), hold on
plot(fcGrid,d1_fd_no(fcGrid),'k')
h0 = errorbar([fc(1:3) fc(5:7)],[prob_d1(1:3) prob_d1(5:7)],[er_d1(1:3) er_d1(5:7)],'o'); set(h0,'Color','k','MarkerFaceColor','k','MarkerSize',1e-10,'LineWidth',2)
plot(fcGrid,far_fd_100(fcGrid),'Color','b')
h1 = errorbar(fc,far_100,erd_100,'o'); set(h1,'Color','b','MarkerFaceColor','b','MarkerSize',1e-10,'LineWidth',2)
plot(fcGrid,far_fd_200(fcGrid),'Color','m')
h2 = errorbar(fc,far_200,erd_200,'o'); set(h2,'Color','m','MarkerFaceColor',[180 0 180]/255,'MarkerSize',1e-10,'LineWidth',2)
plot(fcGrid,far_fd_300(fcGrid),'Color','r')
h3 = errorbar(fc,far_300,erd_300,'o'); set(h3,'Color','r','MarkerFaceColor','r','MarkerSize',1e-10,'LineWidth',2)
errorbar_tick(h0,0); errorbar_tick(h1,0); errorbar_tick(h2,0); errorbar_tick(h3,0); 
axis([100 300 -0.05 1]), set(gca,'YTick',[0 0.5 1])

xlabel('Comparison frequency (Hz)')
ylabel({'P(Comparison', 'judged hihger)'})


% Middle position.
subplot(232), hold on
plot(fcGrid,d1_fd_no(fcGrid),'k')
h0 = errorbar([fc(1:3) fc(5:7)],[prob_d1(1:3) prob_d1(5:7)],[er_d1(1:3) er_d1(5:7)],'o'); set(h0,'Color','k','MarkerFaceColor','k','MarkerSize',1e-10,'LineWidth',2)
plot(fcGrid,middle_fd_100(fcGrid),'Color','b')
h1 = errorbar(fc,middle_100,erm_100,'o'); set(h1,'Color','b','MarkerFaceColor','b','MarkerSize',1e-10,'LineWidth',2)
plot(fcGrid,middle_fd_200(fcGrid),'Color','m')
h2 = errorbar(fc,middle_200,erm_200,'o'); set(h2,'Color','m','MarkerFaceColor',[180 0 180]/255,'MarkerSize',1e-10,'LineWidth',2)
plot(fcGrid,middle_fd_300(fcGrid),'Color','r')
h3 = errorbar(fc,middle_300,erm_300,'o'); set(h3,'Color','r','MarkerFaceColor','r','MarkerSize',1e-10,'LineWidth',2)
errorbar_tick(h0,0); errorbar_tick(h1,0); errorbar_tick(h2,0); errorbar_tick(h3,0); 
axis([100 300 -0.05 1]), set(gca,'YTick',[0 0.5 1])

% Near position.
subplot(233), hold on
plot(fcGrid,d1_fd_no(fcGrid),'k')
h0 = errorbar([fc(1:3) fc(5:7)],[prob_d1(1:3) prob_d1(5:7)],[er_d1(1:3) er_d1(5:7)],'o'); set(h0,'Color','k','MarkerFaceColor','k','MarkerSize',1e-10,'LineWidth',2)
plot(fcGrid,near_fd_100(fcGrid),'Color','b')
h1 = errorbar(fc,near_100,erp_100,'o'); set(h1,'Color','b','MarkerFaceColor','b','MarkerSize',1e-10,'LineWidth',2)
plot(fcGrid,near_fd_200(fcGrid),'Color','m')
h2 = errorbar(fc,near_200,erp_200,'o'); set(h2,'Color','m','MarkerFaceColor',[180 0 180]/255,'MarkerSize',1e-10,'LineWidth',2)
plot(fcGrid,near_fd_300(fcGrid),'Color','r')
h3 = errorbar(fc,near_300,erp_300,'o'); set(h3,'Color','r','MarkerFaceColor','r','MarkerSize',1e-10,'LineWidth',2)
errorbar_tick(h0,0); errorbar_tick(h1,0); errorbar_tick(h2,0); errorbar_tick(h3,0); 

%%
% PSE and JND.
% Fit Gaussian CDF for all stimulus conditions for every subject.
% Extract PSE and JND from the fitted function.
for i = 1:nSubject
    % Fits for D1 and D2 only.
    D1only = fit( [fc(1:3) fc(5:7)]', p_d1(i,[1:3 5:7])', ft, opts );
    D2only = fit( [fc(1:3) fc(5:7)]', p_d2(i,[1:3 5:7])', ft, opts );
    % Fits for near position
    near_fd_100 = fit( fc', p_100(i,:)', ft, opts );
    near_fd_200 = fit( fc', p_200(i,:)', ft, opts );
    near_fd_300 = fit( fc', p_300(i,:)', ft, opts );
    % Fits for middle position.
    middle_fd_100 = fit( fc', m_100(i,:)', ft, opts );
    middle_fd_200 = fit( fc', m_200(i,:)', ft, opts );
    middle_fd_300 = fit( fc', m_300(i,:)', ft, opts );
    % Fits for far position.
    far_fd_100 = fit( fc', d_100(i,:)', ft, opts );
    far_fd_200 = fit( fc', d_200(i,:)', ft, opts );
    far_fd_300 = fit( fc', d_300(i,:)', ft, opts );
    % Extract PSE, indeed the change in PSE with respect to baseline PSE.
    mu_100(i,:) = 1000 * ([far_fd_100.mu middle_fd_100.mu near_fd_100.mu] - D1only.mu);
    mu_200(i,:) = 1000 * ([far_fd_200.mu middle_fd_200.mu near_fd_200.mu] - D1only.mu);
    mu_300(i,:) = 1000 * ([far_fd_300.mu middle_fd_300.mu near_fd_300.mu] - D1only.mu);
    % Extract change in JND with respect to baseline JND.
    sigma_100(i,:) = 1000 * ([far_fd_100.sigma middle_fd_100.sigma near_fd_100.sigma]-D1only.sigma);
    sigma_200(i,:) = 1000 * ([far_fd_200.sigma middle_fd_200.sigma near_fd_200.sigma]-D1only.sigma);
    sigma_300(i,:) = 1000 * ([far_fd_300.sigma middle_fd_300.sigma near_fd_300.sigma]-D1only.sigma);
end

% Combine PSE for each position; also the standard error of mean.
delMU = [mean(mu_100); mean(mu_200); mean(mu_300)];
semMU = [std(mu_100); std(mu_200); std(mu_300)]/sqrt(nSubject);
% Combine JND for each position; also the standard error of mean.
delSIGMA = [mean(sigma_100); mean(sigma_200); mean(sigma_300)];
semSIGMA = [std(sigma_100); std(sigma_200); std(sigma_300)]/sqrt(nSubject);


%% Plot PSE/Bias
% Group by hand locations
subplot(223), X = smartbar((delMU),0.9,1,{'Far','Middle','Near'},c); hold on
% Individual PSe
for ii = 1:3
    plot(4*(ii-1),mu_100(:,ii),'o','MarkerFaceColor','none','MarkerEdgeColor',[0.6 0.6 0.6])
    plot(4*(ii-1)+1,mu_200(:,ii),'o','MarkerFaceColor','none','MarkerEdgeColor',[0.6 0.6 0.6])
    plot(4*(ii-1)+2,mu_300(:,ii),'o','MarkerFaceColor','none','MarkerEdgeColor',[0.6 0.6 0.6])
end

% Error bars
eMU = errorbar(X(:),(delMU(:)),semMU(:),'.'); box off
errorbar_tick(eMU,25)
set(eMU,'Color','k','MarkerFaceColor','k','MarkerSize',1,'LineWidth',1.5)
ylabel('\DeltaPSE (Hz)')
ylim([-65 80])
set(gca,'YTick',[-60 0 80])

text(0,75,'100Hz','FontSize',14,'Color','b', 'HorizontalAlignment','left')
text(0,65,'200Hz','FontSize',14,'Color','m', 'HorizontalAlignment','left')
text(0,55,'300Hz','FontSize',14,'Color','r', 'HorizontalAlignment','left')


%% Plot of JND/Threshold
% Group by hand locations
subplot(224), X = smartbar(delSIGMA,0.9,1,{'Far','Middle','Near'},c); hold on
% Individual JND
for ii = 1:3
    plot(4*(ii-1),sigma_100(:,ii),'o','MarkerFaceColor','none','MarkerEdgeColor',[0.6 0.6 0.6])
    plot(4*(ii-1)+1,sigma_200(:,ii),'o','MarkerFaceColor','none','MarkerEdgeColor',[0.6 0.6 0.6])
    plot(4*(ii-1)+2,sigma_300(:,ii),'o','MarkerFaceColor','none','MarkerEdgeColor',[0.6 0.6 0.6])
end

% Error bars
eSIGMA = errorbar(X(:),delSIGMA(:),semSIGMA(:),'.'); box off
errorbar_tick(eSIGMA,25)
set(eSIGMA,'Color','k','MarkerFaceColor','k','MarkerSize',1,'LineWidth',1.5)
ylabel('\DeltaJND (Hz)')
ylim([-10 55])
set(gca,'YTick',[0 55])
