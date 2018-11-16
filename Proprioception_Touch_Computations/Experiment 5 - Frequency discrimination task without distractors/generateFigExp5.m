clc; clear;

% Set some cool plot/graphics environments!
c =[ 0 0 1;...
     1 0 0;...
     ];
set(0,'defaultAxesColorOrder',c);
set(0,'defaultlinelinewidth',2)
set(0,'DefaultAxesFontSize',16)

subject = {'1','2','3','4','5','6','7','8'};

subjects = strcat('propcut_freq_noDistractor',subject,'.mat');

freq = [100, 140, 180, 220, 260, 300];

n = numel(subjects);

for subj = 1:n
    load(subjects{subj})
    for i= 1:6
        mean_near_cond(subj,i) = mean( data{1,1}{1,i}.prob );
        mean_middle_cond (subj,i) = mean( data{1,2}{1,i}.prob );
        mean_far_cond (subj,i) = mean( data{1,3}{1,i}.prob );
    end
end

% Group mean
groupmean_near = mean(mean_near_cond);
groupmean_middle = mean(mean_middle_cond);
groupmean_far = mean(mean_far_cond);

% Standard error of mean (sem)
groupsem_near = std(mean_near_cond)/sqrt(n);
groupsem_middle = std(mean_middle_cond)/sqrt(n);
groupsem_far = std(mean_far_cond)/sqrt(n);


%Define Gaussian cdf
ft = fittype( '0.5*(1+erf((x-mu)/(sigma*sqrt(2))))', 'independent', 'x');
opts = fitoptions( ft ); opts.Display = 'Off'; opts.Lower = [0 0]; 
opts.StartPoint = [200 50]; opts.Upper = [Inf Inf];

% Group psychometric curves
ft_group_near = fit( freq', groupmean_near', ft, opts );
ft_group_middle = fit( freq', groupmean_middle', ft, opts );
ft_group_far = fit( freq', groupmean_far', ft, opts );

%Plotting smooth curves at group level
finp = 100:5:300;
inp_nearcond = ft_group_near(finp);
inp_middlecond = ft_group_middle(finp);
inp_farcond = ft_group_far(finp);


figure(1), subplot(221)
plot(finp,inp_nearcond,'Color','b'); hold on
plot(finp,inp_middlecond,'Color','g'); hold on
plot(finp,inp_farcond,'Color','r'); hold on

xlabel('Frequency (Hz)');
ylabel('Probability');

% Plot sem on the same figure
eb1 = errorbar(freq,groupmean_near,groupsem_near,'o');
eb2 = errorbar(freq,groupmean_middle,groupsem_middle,'o');
eb3 = errorbar(freq,groupmean_far,groupsem_far,'o');

set(eb1,'color','b','MarkerFaceColor','b','MarkerSize',1e-10); errorbar_tick(eb1,0);
set(eb2,'color','g','MarkerFaceColor','g','MarkerSize',1e-10); errorbar_tick(eb2,0);
set(eb3,'color','r','MarkerFaceColor','r','MarkerSize',1e-10); errorbar_tick(eb3,0);

box off
axis([freq(1) freq(end) 0 1])


% Extracting bias and sensitivity
for subj = 1:n
    ft_nearcond = fit( freq', mean_near_cond(subj,:)', ft, opts );
    bias_nearcond(subj) = ft_nearcond.mu ;
    sensitivity_nearcond(subj) = ft_nearcond.sigma;
    
    ft_middlecond = fit( freq', mean_middle_cond(subj,:)', ft, opts );
    bias_middlecond(subj) = ft_middlecond.mu ;
    sensitivity_middlecond(subj) = ft_middlecond.sigma;
    
    ft_farcond = fit( freq', mean_far_cond(subj,:)', ft, opts );
    bias_farcond(subj) = ft_farcond.mu ;
    sensitivity_farcond(subj) = ft_farcond.sigma;
end

% Mean bias at group level
groupmean_bias_near = mean(bias_nearcond);
groupmean_bias_middle = mean(bias_middlecond);
groupmean_bias_far = mean(bias_farcond);
groupmean_bias = [groupmean_bias_far groupmean_bias_middle groupmean_bias_near];

groupsem_bias_near = std(bias_nearcond)/sqrt(n);
groupsem_bias_middle = std(bias_middlecond)/sqrt(n);
groupsem_bias_far = std(bias_farcond)/sqrt(n);
groupsem_bias = [groupsem_bias_far groupsem_bias_middle groupsem_bias_near];

% Plot of average PSE
position = 1:3;
subplot(2,2,3)
% X = smartbar((groupmean_bias)); hold on
bar(position, groupmean_bias,'FaceColor','none'), hold on
eb = errorbar(position,groupmean_bias, groupsem_bias,'ok');
set(eb,'MarkerFaceColor','k')
errorbar_tick(eb,20);
box off
ylabel('PSE')
set(gca,'XTickLabel',{'Far','Middle','Near'})
axis([0.5 3.5 150 220])

% Mean sensitivity at group level
groupmean_sensitivity_near = mean(sensitivity_nearcond);
groupmean_sensitivity_middle = mean(sensitivity_middlecond);
groupmean_sensitivity_far = mean(sensitivity_farcond);
groupmean_sensitivity = [groupmean_sensitivity_far groupmean_sensitivity_middle groupmean_sensitivity_near];

groupsem_sensitivity_near = std(sensitivity_nearcond)/sqrt(n);
groupsem_sensitivity_middle = std(sensitivity_middlecond)/sqrt(n);
groupsem_sensitivity_far = std(sensitivity_farcond)/sqrt(n);
groupsem_sensitivity = [groupsem_sensitivity_far groupsem_sensitivity_middle groupsem_sensitivity_near];

% Plot of average JND
subplot(2,2,4)
bar(position, groupmean_sensitivity,'FaceColor','none'); hold on
eb = errorbar(position,groupmean_sensitivity, groupsem_sensitivity,'ok');
set(eb,'MarkerFaceColor','k')
errorbar_tick(eb,20);
axis([0.5 3.5 0 80])

% PSE and JND for individual subject
pse = [bias_nearcond; bias_middlecond; bias_farcond];
jnd = [sensitivity_nearcond; sensitivity_middlecond; sensitivity_farcond];

PSE = [pse(3,:); pse(2,:); pse(1,:)];
JND = [jnd(3,:); jnd(2,:); jnd(1,:)];

for i = 1:size(jnd,2)
    subplot(2,2,3)
    plot(position,PSE(:,i),'ok','LineWidth',0.5)
    subplot(2,2,4)
    plot(position,JND(:,i),'ok','LineWidth',0.5)
end

box off
ylabel('JND')
set(gca,'XTickLabel',{'Far','Middle','Near'})








