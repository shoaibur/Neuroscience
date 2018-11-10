clc; clear;

% Read data as a structure
load('freqFullData.mat')

% Values used in the experiment
standard = 200;
% comparison frequency
fc = repmat([100 140 180 200 220 260 300], 9, 1);
% hand positions
prop = repmat([0*ones(3,1); 0.5*ones(3,1); ones(3,1)], 1, 7);
% distractor frequency
fd = repmat( repmat([100;200;300],3,7), 1, 1);
% number of repeats
n = repmat([20 22 24 24 24 22 20], 9, 1);

% Initialization for parameter search
c_init_si = [
    0.4524      -0.8804      0.49372      0.95277        18097
    0.4695       -0.763      0.50442       1.0561        18446
    0.4793      -1.2629      0.53364      0.56515        13806
    0.5         -1.625       0.52904      0.52594        14026
    0.4539      -0.7704      0.48656       1.0659        22936
    0.4541      -1.0743      0.52087      0.69247        17071
    0.4671      -1.3744      0.45694      0.59343        18429
    0.4693      -0.7823      0.55619       1.0338        17290];
%
% 2D frequency- and location-based modulation function
afs = @(x1,x2,x3,dp,df) x1 * exp(-( dp/x2 + abs(df).^2/x3 ));
afsn = @(x1,x2,dp,df) exp(-( dp/x1 + abs(df).^2/x2 ));
% 1D frequency-based modulation function
afsf = @(x1,x2,df) x1 * exp(-( abs(df).^2/x2 ));
afsfn = @(x1,df) exp(-( abs(df).^2/x1 ));
% 1D location-based modulation function
afsp = @(x1,x2,dp) x1 * exp(-( dp/x2 ));
afspn = @(x1,dp) exp(-( dp/x1 ));

for subj = 1:8
    % Single subject    
    pseD1 = out.pseD1(9*subj-8);
    jndD1 = out.jndD1(9*subj-8);
    jndD2 = out.jndD2(9*subj-8);
    subj_cp = out.choiceprob(9*subj-8:9*subj,:);

    % Count of subject choice probability
    X = subj_cp .* n;
    
    % Parameter free model
    llfree(subj) = sum(sum(...
        X.*log(0.5 * ( 1 + erf(...
        (fc - repmat(out.pse(9*subj-8:9*subj),1,7) ) ./...
        sqrt( 2 * (repmat(out.jnd(9*subj-8:9*subj),1,7)).^2 )...
        ) ) ) + ...
        (n-X).*log(0.5 * ( 1 - erf(...
        (fc - repmat(out.pse(9*subj-8:9*subj),1,7) ) ./...
        sqrt( 2 * (repmat(out.jnd(9*subj-8:9*subj),1,7)).^2 )...
        ) ) )...
        ));
    
    % FP+FP, without normalization
    % 2D-2D-aa
    wT = @(c,prop,df) c(4);
    wD = @(c,prop,df) afs(c(4),c(5),c(6),prop,df);
    f2D2Daa = @(c,fc,prop,standard,fd,jndD1,jndD2)...
        -sum(sum(...
        X.*log(0.5 * ( 1 + erf(...
        (fc - wT(c,prop,(standard-fd))*standard - wD(c,prop,(standard-fd)).*fd ) ./...
        sqrt( 2 * ( (wT(c,prop,(standard-fd))*jndD1/c(1)).^2 + (wD(c,prop,(standard-fd)).*jndD2./afs(c(1),c(2),c(3),prop,(standard-fd))).^2 ) )...
        ) ) ) + ...
        (n-X).*log(0.5 * ( 1 - erf(...
        (fc - wT(c,prop,(standard-fd))*standard - wD(c,prop,(standard-fd)).*fd ) ./...
        sqrt( 2 * ( (wT(c,prop,(standard-fd))*jndD1/c(1)).^2 + (wD(c,prop,(standard-fd)).*jndD2./afs(c(1),c(2),c(3),prop,(standard-fd))).^2 ) )...
        ) ) )...
        ));
    
    fun = @(c)f2D2Daa(c,fc,prop,standard,fd,jndD1,jndD2);
    c2D2Daa(subj,:) = fminsearch( fun, [c_init_si(subj,3) 1000 c_init_si(subj,5) c_init_si(subj,3) c_init_si(subj,4) 1000000], optimset('MaxFunEvals',100) );
    ll_2D2Daa(subj) = -f2D2Daa(c2D2Daa(subj,:),fc,prop,standard,fd,jndD1,jndD2);
    %}
    
    
    % FP+FP, with normalization
    % 2D-2D-ra
    wT = @(c,prop,df) 1 ./ (c(4) + afsn(c(5),c(6),prop,df));
    wD = @(c,prop,df) afsn(c(5),c(6),prop,df) ./ (c(4) + afsn(c(5),c(6),prop,df));
    f2D2Dra = @(c,fc,prop,standard,fd,jndD1,jndD2)...
        -sum(sum(...
        X.*log(0.5 * ( 1 + erf(...
        (fc - wT(c,prop,(standard-fd))*standard - wD(c,prop,(standard-fd)).*fd ) ./...
        sqrt( 2 * ( (wT(c,prop,(standard-fd)).*jndD1/c(1)).^2 + (wD(c,prop,(standard-fd)).*jndD2./afs(c(1),c(2),c(3),prop,(standard-fd))).^2 ) )...
        ) ) ) + ...
        (n-X).*log(0.5 * ( 1 - erf(...
        (fc - wT(c,prop,(standard-fd))*standard - wD(c,prop,(standard-fd)).*fd ) ./...
        sqrt( 2 * ( (wT(c,prop,(standard-fd)).*jndD1/c(1)).^2 + (wD(c,prop,(standard-fd)).*jndD2./afs(c(1),c(2),c(3),prop,(standard-fd))).^2 ) )...
        ) ) )...
        ));
    
    fun = @(c)f2D2Dra(c,fc,prop,standard,fd,jndD1,jndD2);
    c2D2Dra(subj,:) = fminsearch( fun, [c_init_si(subj,[3 4 5]) 1 c_init_si(subj,[4 5])], optimset('MaxFunEvals',100) );
    ll_2D2Dra(subj) = -f2D2Dra(c2D2Dra(subj,:),fc,prop,standard,fd,jndD1,jndD2);
    %}

    
    % FP:P, without normalization
    % 2D:1Dp-aa
    wT = @(c,prop,df) c(4);
    wD = @(c,prop,df) afsp(c(4),c(3),prop);
    f2D1Dpaa = @(c,fc,prop,standard,fd,jndD1,jndD2)...
        -sum(sum(...
        X.*log(0.5 * ( 1 + erf(...
        (fc - wT(c,prop,(standard-fd))*standard - wD(c,prop,(standard-fd)).*fd ) ./...
        sqrt( 2 * ( (wT(c,prop,(standard-fd))*jndD1/c(1)).^2 + (wD(c,prop,(standard-fd)).*jndD2./afs(c(1),c(2),c(3),prop,(standard-fd))).^2 ) )...
        ) ) ) + ...
        (n-X).*log(0.5 * ( 1 - erf(...
        (fc - wT(c,prop,(standard-fd))*standard - wD(c,prop,(standard-fd)).*fd ) ./...
        sqrt( 2 * ( (wT(c,prop,(standard-fd))*jndD1/c(1)).^2 + (wD(c,prop,(standard-fd)).*jndD2./afs(c(1),c(2),c(3),prop,(standard-fd))).^2 ) )...
        ) ) )...
        ));
    
    fun = @(c)f2D1Dpaa(c,fc,prop,standard,fd,jndD1,jndD2);
    c2D1Dpaa(subj,:) = fminsearch( fun, [c_init_si(subj,[3 4 5]) c_init_si(subj,3)], optimset('MaxFunEvals',100) );
    ll_2D1Dpaa(subj) = -f2D1Dpaa(c2D1Dpaa(subj,:),fc,prop,standard,fd,jndD1,jndD2);
    %}
    
    
    % FP:P, with normalization
    % 2D:1Dp-ra
    wT = @(c,prop,df) 1 ./ (c(4) + afspn(c(3),prop));
    wD = @(c,prop,df) afspn(c(3),prop) ./ (c(4) + afspn(c(3),prop));
    f2D1Dpra = @(c,fc,prop,standard,fd,jndD1,jndD2)...
        -sum(sum(...
        X.*log(0.5 * ( 1 + erf(...
        (fc - wT(c,prop,(standard-fd))*standard - wD(c,prop,(standard-fd)).*fd ) ./...
        sqrt( 2 * ( (wT(c,prop,(standard-fd)).*jndD1/c(1)).^2 + (wD(c,prop,(standard-fd)).*jndD2./afs(c(1),c(2),c(3),prop,(standard-fd))).^2 ) )...
        ) ) ) + ...
        (n-X).*log(0.5 * ( 1 - erf(...
        (fc - wT(c,prop,(standard-fd))*standard - wD(c,prop,(standard-fd)).*fd ) ./...
        sqrt( 2 * ( (wT(c,prop,(standard-fd)).*jndD1/c(1)).^2 + (wD(c,prop,(standard-fd)).*jndD2./afs(c(1),c(2),c(3),prop,(standard-fd))).^2 ) )...
        ) ) )...
        ));
    
    fun = @(c)f2D1Dpra(c,fc,prop,standard,fd,jndD1,jndD2);
    c2D1Dpra(subj,:) = fminsearch( fun, [c_init_si(subj,[3 4 5]) 1], optimset('MaxFunEvals',100) );
    ll_2D1Dpra(subj) = -f2D1Dpra(c2D1Dpra(subj,:),fc,prop,standard,fd,jndD1,jndD2);
    %}
    
    
    % FP+P, without normalization
    % 2D-1Dp-aa
    wT = @(c,prop,df) c(4);
    wD = @(c,prop,df) afsp(c(4),c(5),prop);
    f2D1Daa = @(c,fc,prop,standard,fd,jndD1,jndD2)...
        -sum(sum(...
        X.*log(0.5 * ( 1 + erf(...
        (fc - wT(c,prop,(standard-fd))*standard - wD(c,prop,(standard-fd)).*fd ) ./...
        sqrt( 2 * ( (wT(c,prop,(standard-fd))*jndD1/c(1)).^2 + (wD(c,prop,(standard-fd)).*jndD2./afs(c(1),c(2),c(3),prop,(standard-fd))).^2 ) )...
        ) ) ) + ...
        (n-X).*log(0.5 * ( 1 - erf(...
        (fc - wT(c,prop,(standard-fd))*standard - wD(c,prop,(standard-fd)).*fd ) ./...
        sqrt( 2 * ( (wT(c,prop,(standard-fd))*jndD1/c(1)).^2 + (wD(c,prop,(standard-fd)).*jndD2./afs(c(1),c(2),c(3),prop,(standard-fd))).^2 ) )...
        ) ) )...
        ));
    
    fun = @(c)f2D1Daa(c,fc,prop,standard,fd,jndD1,jndD2);
    c2D1Daa(subj,:) = fminsearch( fun, [c_init_si(subj,[3 4 5]) c_init_si(subj,[3 4])], optimset('MaxFunEvals',100) );
    ll_2D1Daa(subj) = -f2D1Daa(c2D1Daa(subj,:),fc,prop,standard,fd,jndD1,jndD2);
    %}
    
    
    % FP+P, with normalization
    % 2D-1Dp-ra
    wT = @(c,prop,df) 1 ./ (c(4) + afspn(c(5),prop));
    wD = @(c,prop,df) afspn(c(5),prop) ./ (c(4) + afspn(c(5),prop));
    f2D1Dra = @(c,fc,prop,standard,fd,jndD1,jndD2)...
        -sum(sum(...
        X.*log(0.5 * ( 1 + erf(...
        (fc - wT(c,prop,(standard-fd))*standard - wD(c,prop,(standard-fd)).*fd ) ./...
        sqrt( 2 * ( (wT(c,prop,(standard-fd)).*jndD1/c(1)).^2 + (wD(c,prop,(standard-fd)).*jndD2./afs(c(1),c(2),c(3),prop,(standard-fd))).^2 ) )...
        ) ) ) + ...
        (n-X).*log(0.5 * ( 1 - erf(...
        (fc - wT(c,prop,(standard-fd))*standard - wD(c,prop,(standard-fd)).*fd ) ./...
        sqrt( 2 * ( (wT(c,prop,(standard-fd)).*jndD1/c(1)).^2 + (wD(c,prop,(standard-fd)).*jndD2./afs(c(1),c(2),c(3),prop,(standard-fd))).^2 ) )...
        ) ) )...
        ));
    
    fun = @(c)f2D1Dra(c,fc,prop,standard,fd,jndD1,jndD2);
    c2D1Dra(subj,:) = fminsearch( fun, [c_init_si(subj,[3 4 5]) 1 c_init_si(subj,4)], optimset('MaxFunEvals',100) );
    ll_2D1Dra(subj) = -f2D1Dra(c2D1Dra(subj,:),fc,prop,standard,fd,jndD1,jndD2);
    %}
    
    
    
    % FP, optimal cue integration
    % 2D-rr
    wT = @(c,prop,df,jndD1,jndD2) (jndD2./afs(c(1),c(2),c(3),prop,df)).^2 ./ ( (jndD1/c(1)).^2 + (jndD2./afs(c(1),c(2),c(3),prop,df)).^2 );
    wD = @(c,prop,df,jndD1,jndD2) (jndD1/c(1)).^2 ./ ( (jndD1/c(1)).^2 + (jndD2./afs(c(1),c(2),c(3),prop,df)).^2 );
    f2Drr = @(c,fc,prop,standard,fd,jndD1,jndD2)...
        -sum(sum(...
        X.*log(0.5 * ( 1 + erf(...
        (fc - wT(c,prop,(standard-fd),jndD1,jndD2)*standard - wD(c,prop,(standard-fd),jndD1,jndD2).*fd ) ./...
        sqrt( 2 * ( wT(c,prop,(standard-fd),jndD1,jndD2).^2 .* (jndD2./afs(c(1),c(2),c(3),prop,(standard-fd))).^2 + wD(c,prop,(standard-fd),jndD1,jndD2).^2 .* (jndD1/c(1)).^2 ) )...
        ) ) ) + ...
        (n-X).*log(0.5 * ( 1 - erf(...
        (fc - wT(c,prop,(standard-fd),jndD1,jndD2)*standard - wD(c,prop,(standard-fd),jndD1,jndD2).*fd ) ./...
        sqrt( 2 * ( wT(c,prop,(standard-fd),jndD1,jndD2).^2 .* (jndD2./afs(c(1),c(2),c(3),prop,(standard-fd))).^2 + wD(c,prop,(standard-fd),jndD1,jndD2).^2 .* (jndD1/c(1)).^2 ) )...
        ) ) )...
        ));
    
    fun = @(c)f2Drr(c,fc,prop,standard,fd,jndD1,jndD2);
    c2Drr(subj,:) = fminsearch( fun, c_init_si(subj,3:5), optimset('MaxFunEvals',100) );
    ll_2Drr(subj) = -f2Drr(c2Drr(subj,:),fc,prop,standard,fd,jndD1,jndD2);
    %}
    
    
    % FP:FP, without normalization
    % 2D-aa
    wT = @(c,prop,df) c(1);
    wD = @(c,prop,df) afs(c(1),c(2),c(3),prop,df);
    f2Daa = @(c,fc,prop,standard,fd,jndD1,jndD2)...
        -sum(sum(...
        X.*log(0.5 * ( 1 + erf(...
        (fc - wT(c,prop,(standard-fd))*standard - wD(c,prop,(standard-fd)).*fd ) ./...
        sqrt( 2 * ( jndD1.^2 + jndD2.^2 ) )...
        ) ) ) + ...
        (n-X).*log(0.5 * ( 1 - erf(...
        (fc - wT(c,prop,(standard-fd))*standard - wD(c,prop,(standard-fd)).*fd ) ./...
        sqrt( 2 * ( jndD1.^2 + jndD2.^2 ) )...
        ) ) )...
        ));
    
    fun = @(c)f2Daa(c,fc,prop,standard,fd,jndD1,jndD2);
    c2Daa(subj,:) = fminsearch( fun, c_init_si(subj,3:5), optimset('MaxFunEvals',100) );
    ll_2Daa(subj) = -f2Daa(c2Daa(subj,:),fc,prop,standard,fd,jndD1,jndD2);
    %}
    
    
    
    % FP:FP, with normalization
    % 2D-ra
    wT = @(c,prop,df) 1./(c(4)+afsn(c(2),c(3),prop,df));
    wD = @(c,prop,df) afsn(c(2),c(3),prop,df)./(c(4)+afsn(c(2),c(3),prop,df));
    f2Dra = @(c,fc,prop,standard,fd,jndD1,jndD2)...
        -sum(sum(...
        X.*log(0.5 * ( 1 + erf(...
        (fc - wT(c,prop,(standard-fd))*standard - wD(c,prop,(standard-fd)).*fd ) ./...
        sqrt( 2 * ( ( wT(c,prop,(standard-fd)) .* jndD1 / c(1) ).^2 + ( wD(c,prop,(standard-fd)) .* jndD2./afs(c(1),c(2),c(3),prop,(standard-fd)) ).^2 ) )...
        ) ) ) + ...
        (n-X).*log(0.5 * ( 1 - erf(...
        (fc - wT(c,prop,(standard-fd))*standard - wD(c,prop,(standard-fd)).*fd ) ./...
        sqrt( 2 * ( ( wT(c,prop,(standard-fd)) .* jndD1 / c(1) ).^2 + ( wD(c,prop,(standard-fd)) .* jndD2./afs(c(1),c(2),c(3),prop,(standard-fd)) ).^2 ) )...
        ) ) )...
        ));
    
    fun = @(c)f2Dra(c,fc,prop,standard,fd,jndD1,jndD2);
    c2Dra(subj,:) = fminsearch( fun, [c_init_si(subj,3:5) 1], optimset('MaxFunEvals',100) );
    ll_2Dra(subj) = -f2Dra(c2Dra(subj,:),fc,prop,standard,fd,jndD1,jndD2);
    %}
    
    
    
    % P, optimal cue integration
    % 1Dp-rr
    wT = @(c,dp,jndD1,jndD2) (jndD2./afsp(c(1),c(2),dp)).^2 ./ ( (jndD1/c(1)).^2 + (jndD2./afsp(c(1),c(2),dp)).^2);
    wD = @(c,dp,jndD1,jndD2) (jndD1/c(1)).^2 ./ ( (jndD1/c(1)).^2 + (jndD2./afsp(c(1),c(2),dp)).^2);
    f1Dprr = @(c,fc,prop,standard,fd,jndD1,jndD2)...
        -sum(sum(...
        X.*log(0.5 * ( 1 + erf(...
        (fc - wT(c,prop,jndD1,jndD2)*standard - wD(c,prop,jndD1,jndD2).*fd ) ./...
        sqrt( 2 * ( (wT(c,prop,jndD1,jndD2).*jndD1/c(1)).^2 + (wD(c,prop,jndD1,jndD2).*jndD2./afsp(c(1),c(2),prop)).^2) )...
        ) ) ) + ...
        (n-X).*log(0.5 * ( 1 - erf(...
        (fc - wT(c,prop,jndD1,jndD2)*standard - wD(c,prop,jndD1,jndD2).*fd ) ./...
        sqrt( 2 * ( (wT(c,prop,jndD1,jndD2).*jndD1/c(1)).^2 + (wD(c,prop,jndD1,jndD2).*jndD2./afsp(c(1),c(2),prop)).^2) )...
        ) ) )...
        ));
    
    fun = @(c)f1Dprr(c,fc,prop,standard,fd,jndD1,jndD2);
    c1Dprr(subj,:) = fminsearch( fun, c_init_si(subj,3:4), optimset('MaxFunEvals',100) );
    ll_1Dprr(subj) = -f1Dprr(c1Dprr(subj,:),fc,prop,standard,fd,jndD1,jndD2);
    %}
    
    
    % P, without normalization
    % 1Dp-aa
    f1Dpaa = @(c,fc,prop,standard,fd,jndD1,jndD2)...
        -sum(sum(...
        X.*log(0.5 * ( 1 + erf(...
        (fc - ( c(1) )*standard - ( afsp(c(1),c(2),prop) ).*fd ) ./...
        sqrt( 2 * ( jndD1.^2 + jndD2.^2 ) )...
        ) ) ) + ...
        (n-X).*log(0.5 * ( 1 - erf(...
        (fc - ( c(1) )*standard - ( afsp(c(1),c(2),prop) ).*fd ) ./...
        sqrt( 2 * ( jndD1.^2 + jndD2.^2 ) )...
        ) ) )...
        ));
    
    fun = @(c)f1Dpaa(c,fc,prop,standard,fd,jndD1,jndD2);
    c1Dpaa(subj,:) = fminsearch( fun, c_init_si(subj,3:4), optimset('MaxFunEvals',100) );
    ll_1Dpaa(subj) = -f1Dpaa(c1Dpaa(subj,:),fc,prop,standard,fd,jndD1,jndD2);
    %}
    

    % P, with normalization
    % 1Dp-ra
    wT = @(c,prop) c(1)./(c(3)+c(1)+afsp(c(1),c(2),prop));
    wD = @(c,prop) afsp(c(1),c(2),prop)./(c(3)+c(1)+afsp(c(1),c(2),prop));
    f1Dpra = @(c,fc,prop,standard,fd,jndD1,jndD2)...
        -sum(sum(...
        X.*log(0.5 * ( 1 + erf(...
        (fc - wT(c,prop)*standard - wD(c,prop).*fd ) ./...
        sqrt( 2 * ( (wT(c,prop).*jndD1/c(1)).^2 + (wD(c,prop).*jndD2./afsp(c(1),c(2),prop)).^2) )...
        ) ) ) + ...
        (n-X).*log(0.5 * ( 1 - erf(...
        (fc - wT(c,prop)*standard - wD(c,prop).*fd ) ./...
        sqrt( 2 * ( (wT(c,prop).*jndD1/c(1)).^2 + (wD(c,prop).*jndD2./afsp(c(1),c(2),prop)).^2) )...
        ) ) )...
        ));
    
    fun = @(c)f1Dpra(c,fc,prop,standard,fd,jndD1,jndD2);
    c1Dpra(subj,:) = fminsearch( fun, [c_init_si(subj,3:4) 0], optimset('MaxFunEvals',100) );
    ll_1Dpra(subj) = -f1Dpra(c1Dpra(subj,:),fc,prop,standard,fd,jndD1,jndD2);
    %}
    

    
    % F, optimal cue integration
    % 1Df-rr
    wT = @(c,df,jndD1,jndD2) (jndD2./afsf(c(1),c(2),df)).^2 ./ ( (jndD1/c(1)).^2 + (jndD2./afsf(c(1),c(2),df)).^2);
    wD = @(c,df,jndD1,jndD2) (jndD1/c(1)).^2 ./ ( (jndD1/c(1)).^2 + (jndD2./afsf(c(1),c(2),df)).^2);
    f1Dfrr = @(c,fc,prop,standard,fd,jndD1,jndD2)...
        -sum(sum(...
        X.*log(0.5 * ( 1 + erf(...
        (fc - wT(c,(standard-fd),jndD1,jndD2)*standard - wD(c,(standard-fd),jndD1,jndD2).*fd ) ./...
        sqrt( 2 * ( (wT(c,(standard-fd),jndD1,jndD2).*jndD1/c(1)).^2 + (wD(c,(standard-fd),jndD1,jndD2).*jndD2./afsp(c(1),c(2),(standard-fd))).^2) )...
        ) ) ) + ...
        (n-X).*log(0.5 * ( 1 - erf(...
        (fc - wT(c,(standard-fd),jndD1,jndD2)*standard - wD(c,(standard-fd),jndD1,jndD2).*fd ) ./...
        sqrt( 2 * ( (wT(c,(standard-fd),jndD1,jndD2).*jndD1/c(1)).^2 + (wD(c,(standard-fd),jndD1,jndD2).*jndD2./afsp(c(1),c(2),(standard-fd))).^2) )...
        ) ) )...
        ));
    
    fun = @(c)f1Dfrr(c,fc,prop,standard,fd,jndD1,jndD2);
    c1Dfrr(subj,:) = fminsearch( fun, c_init_si(subj,[3 5]), optimset('MaxFunEvals',100) );
    ll_1Dfrr(subj) = -f1Dfrr(c1Dfrr(subj,:),fc,prop,standard,fd,jndD1,jndD2);
    %}
    
    
    
    % F, without normalization
    % 1Df-aa
    f1Dfaa = @(c,fc,prop,standard,fd,jndD1,jndD2)...
        -sum(sum(...
        X.*log(0.5 * ( 1 + erf(...
        (fc - ( c(1) )*standard - ( afsf(c(1),c(2),(standard-fd)) ).*fd ) ./...
        sqrt( 2 * ( jndD1.^2 + jndD2.^2 ) )...
        ) ) ) + ...
        (n-X).*log(0.5 * ( 1 - erf(...
        (fc - ( c(1) )*standard - ( afsf(c(1),c(2),(standard-fd)) ).*fd ) ./...
        sqrt( 2 * ( jndD1.^2 + jndD2.^2 ) )...
        ) ) )...
        ));
    
    fun = @(c)f1Dfaa(c,fc,prop,standard,fd,jndD1,jndD2);
    c1Dfaa(subj,:) = fminsearch( fun, c_init_si(subj,[3 5]), optimset('MaxFunEvals',100) );
    ll_1Dfaa(subj) = -f1Dfaa(c1Dfaa(subj,:),fc,prop,standard,fd,jndD1,jndD2);
    %}
    
    

    
    % F, with normalization
    % 1Df-ra
    wT = @(c,df) c(1)./(c(3)+c(1)+afsf(c(1),c(2),df));
    wD = @(c,df) afsf(c(1),c(2),df)./(c(3)+c(1)+afsf(c(1),c(2),df));
    f1Dfra = @(c,fc,prop,standard,fd,jndD1,jndD2)...
        -sum(sum(...
        X.*log(0.5 * ( 1 + erf(...
        (fc - wT(c,(standard-fd))*standard - wD(c,(standard-fd)).*fd ) ./...
        sqrt( 2 * ( (wT(c,(standard-fd)).*jndD1/c(1)).^2 + (wD(c,(standard-fd)).*jndD2./afsf(c(1),c(2),(standard-fd))).^2) )...
        ) ) ) + ...
        (n-X).*log(0.5 * ( 1 - erf(...
        (fc - wT(c,(standard-fd))*standard - wD(c,(standard-fd)).*fd ) ./...
        sqrt( 2 * ( (wT(c,(standard-fd)).*jndD1/c(1)).^2 + (wD(c,(standard-fd)).*jndD2./afsf(c(1),c(2),(standard-fd))).^2) )...
        ) ) )...
        ));
    
    fun = @(c)f1Dfra(c,fc,prop,standard,fd,jndD1,jndD2);
    c1Dfra(subj,:) = fminsearch( fun, [c_init_si(subj,[3 5]) 0], optimset('MaxFunEvals',100) );
    ll_1Dfra(subj) = -f1Dfra(c1Dfra(subj,:),fc,prop,standard,fd,jndD1,jndD2);
    %}
    
    
    
    % F+P, without normalization
    % 1Df, 1Dp - aa
    wT = @(c,prop) c(3);
    wD = @(c,prop) afsp(c(3),c(4),prop);
    f1Df1Dpaa = @(c,fc,prop,standard,fd,jndD1,jndD2)...
        -sum(sum(...
        X.*log(0.5 * ( 1 + erf(...
        (fc - wT(c,prop)*standard - wD(c,prop).*fd ) ./...
        sqrt( 2 * ( (wT(c,prop).*jndD1/c(1)).^2 + (wD(c,prop).*jndD2./afsf(c(1),c(2),(standard-fd))).^2 ) )...
        ) ) ) + ...
        (n-X).*log(0.5 * ( 1 - erf(...
        (fc - wT(c,prop)*standard - wD(c,prop).*fd ) ./...
        sqrt( 2 * ( (wT(c,prop).*jndD1/c(1)).^2 + (wD(c,prop).*jndD2./afsf(c(1),c(2),(standard-fd))).^2 ) )...
        ) ) )...
        ));
    
    fun = @(c)f1Df1Dpaa(c,fc,prop,standard,fd,jndD1,jndD2);
    c1Df1Dpaa(subj,:) = fminsearch( fun, c_init_si(subj,[3 5 3 4]), optimset('MaxFunEvals',100) );
    ll_DfDpaa(subj) = -f1Df1Dpaa(c1Df1Dpaa(subj,:),fc,prop,standard,fd,jndD1,jndD2);
    %}
    
    
    
    
    % F+P, with normalization
    % 1Df, 1Dp - ra
    wT = @(c,prop) 1 ./ ( c(4)+afspn(c(3),prop) );
    wD = @(c,prop) afspn(c(3),prop) ./ ( c(4)+afspn(c(3),prop) );
    %wT = @(c,df,jndD1,jndD2) 1 ./ ( (sqrt(jndD1.^2+jndD2.^2)/(c(1)*jndD1))+afsfn(c(3),df) );
    %wD = @(c,df,jndD1,jndD2) afsfn(c(3),df) ./ ( (sqrt(jndD1.^2+jndD2.^2)/(c(1)*jndD1))+afsfn(c(3),df) );
    f1Df1Dpra = @(c,fc,prop,standard,fd,jndD1,jndD2)...
        -sum(sum(...
        X.*log(0.5 * ( 1 + erf(...
        (fc - wT(c,prop)*standard - wD(c,prop).*fd ) ./...
        sqrt( 2 * ( (wT(c,prop).*jndD1/c(1)).^2 + (wD(c,prop).*jndD2./afsf(c(1),c(2),(standard-fd))).^2 ) )...
        ) ) ) + ...
        (n-X).*log(0.5 * ( 1 - erf(...
        (fc - wT(c,prop)*standard - wD(c,prop).*fd ) ./...
        sqrt( 2 * ( (wT(c,prop).*jndD1/c(1)).^2 + (wD(c,prop).*jndD2./afsf(c(1),c(2),(standard-fd))).^2 ) )...
        ) ) )...
        ));
    
    fun = @(c)f1Df1Dpra(c,fc,prop,standard,fd,jndD1,jndD2);
    c1Df1Dpra(subj,:) = fminsearch( fun, [c_init_si(subj,[3 5 4]) 1], optimset('MaxFunEvals',100) );
    ll_DfDpra(subj) = -f1Df1Dpra(c1Df1Dpra(subj,:),fc,prop,standard,fd,jndD1,jndD2);
    %}
    
    
    
    % P+F, without normalization
    % 1Dp, 1Df - aa
    wT = @(c,df) c(3);
    wD = @(c,df) afsf(c(3),c(4),df);
    f1Dp1Dfaa = @(c,fc,prop,standard,fd,jndD1,jndD2)...
        -sum(sum(...
        X.*log(0.5 * ( 1 + erf(...
        (fc - wT(c,(standard-fd))*standard - wD(c,(standard-fd)).*fd ) ./...
        sqrt( 2 * ( (wT(c,(standard-fd)).*jndD1/c(1)).^2 + (wD(c,(standard-fd)).*jndD2./afsp(c(1),c(2),prop)).^2 ) )...
        ) ) ) + ...
        (n-X).*log(0.5 * ( 1 - erf(...
        (fc - wT(c,(standard-fd))*standard - wD(c,(standard-fd)).*fd ) ./...
        sqrt( 2 * ( (wT(c,(standard-fd)).*jndD1/c(1)).^2 + (wD(c,(standard-fd)).*jndD2./afsp(c(1),c(2),prop)).^2 ) )...
        ) ) )...
        ));
    
    fun = @(c)f1Dp1Dfaa(c,fc,prop,standard,fd,jndD1,jndD2);
    c1Dp1Dfaa(subj,:) = fminsearch( fun, c_init_si(subj,[3 4 3 5]), optimset('MaxFunEvals',100) );
    ll_DpDfaa(subj) = -f1Dp1Dfaa(c1Df1Dpaa(subj,:),fc,prop,standard,fd,jndD1,jndD2);
    %}
    
    
    
    % P+F, with normalization
    % 1Dp, 1Df - ra
    wT = @(c,df) 1 ./ ( c(4)+afsfn(c(3),df) );
    wD = @(c,df) afsfn(c(3),df) ./ ( c(4)+afsfn(c(3),df) );
    f1Dp1Dfra = @(c,fc,prop,standard,fd,jndD1,jndD2)...
        -sum(sum(...
        X.*log(0.5 * ( 1 + erf(...
        (fc - wT(c,(standard-fd))*standard - wD(c,(standard-fd)).*fd ) ./...
        sqrt( 2 * ( (wT(c,(standard-fd)).*jndD1/c(1)).^2 + (wD(c,(standard-fd)).*jndD2./afsp(c(1),c(2),prop)).^2 ) )...
        ) ) ) + ...
        (n-X).*log(0.5 * ( 1 - erf(...
        (fc - wT(c,(standard-fd))*standard - wD(c,(standard-fd)).*fd ) ./...
        sqrt( 2 * ( (wT(c,(standard-fd)).*jndD1/c(1)).^2 + (wD(c,(standard-fd)).*jndD2./afsp(c(1),c(2),prop)).^2 ) )...
        ) ) )...
        ));
    
    fun = @(c)f1Dp1Dfra(c,fc,prop,standard,fd,jndD1,jndD2);
    c1Dp1Dfra(subj,:) = fminsearch( fun, [c_init_si(subj,[3 4 5]) 1], optimset('MaxFunEvals',100) );
    ll_DpDrra(subj) = -f1Dp1Dfra(c1Dp1Dfra(subj,:),fc,prop,standard,fd,jndD1,jndD2);
    %}
end

%
clc;
% Aggregate log-likelihood
ll = [ll_2D2Daa; ll_2D2Dra; ll_2D1Daa; ll_2D1Dra; ll_2D1Dpaa; ll_2D1Dpra; ll_2Drr; ll_2Daa; ll_2Dra; ll_1Dprr; ll_1Dpaa; ll_1Dpra; ll_1Dfrr; ll_1Dfaa; ll_1Dfra; ll_DfDpaa; ll_DfDpra; ll_DpDfaa; ll_DpDrra]';
% Compute log-likelihood ratio test
lrt = -2*log(ll ./ repmat(ll(:,17),1,19));
% Number of parameters in each model; need for AIC & BIC computation
numFreeParam = [6 6 5 5, 4 4, 3 3 4, 2 2 3, 2 2 3, 4 4, 4 4];
% May be used for corrected AIC; not used here!
AICcorrection = (2*numFreeParam.^2 + 2*numFreeParam) ./ (9*7-numFreeParam-1);
% Compute AIC and BIC
for subj = 1:8
    [a,b] = aicbic(ll(subj,:),numFreeParam,9*7*ones(1,numel(numFreeParam)));
    aic(subj,:) = a; bic(subj,:) = b;
end


% Compute model probabilities based on estimated c parameters
fc = [100 140 180 200 220 260 300];
prop = [0 0.5 1];
fd = [100 200 300];
row = 0;
for subj = 1:8
    pseD1 = out.pseD1(9*subj - 8);
    jndD1 = out.jndD1(9*subj - 8);
    jndD2 = out.jndD2(9*subj - 8);
    for prop = [0 0.5 1]
        for fd = [100 200 300]
            row = row + 1;
            
            wT = c2D2Daa(subj,4);
            wD = afs(c2D2Daa(subj,4),c2D2Daa(subj,5),c2D2Daa(subj,6),prop,(standard-fd));
            pse_2D2Daa(row) = wT*standard + wD.*fd;
            jnd_2D2Daa(row) = sqrt( (wT*jndD1/c2D2Daa(subj,1)).^2 + (wD.*jndD2./afs(c2D2Daa(subj,1),c2D2Daa(subj,2),c2D2Daa(subj,3),prop,(standard-fd))).^2 );
            mcp_2D2Daa(row,:) = 0.5 * ( 1+erf( (fc - pse_2D2Daa(row) ) ./ sqrt( 2 * jnd_2D2Daa(row)^2 ) ) );
            
            rsscp_2D2Daa(row) = sum( (mcp_2D2Daa(row,:)-out.choiceprob(row,:)).^2 );
            rsspse_2D2Daa(row) = sum( (pse_2D2Daa(row)-out.pse(row)).^2 );
            rssjnd_2D2Daa(row) = sum( (jnd_2D2Daa(row)-out.jnd(row)).^2 );
            
            
            wT = 1 ./ (c2D2Dra(subj,4) + afsn(c2D2Dra(subj,5),c2D2Dra(subj,6),prop,(standard-fd)));
            wD = afsn(c2D2Dra(subj,5),c2D2Dra(subj,6),prop,(standard-fd)) ./ (c2D2Dra(subj,4) + afsn(c2D2Dra(subj,5),c2D2Dra(subj,6),prop,(standard-fd)));
            pse_2D2Dra(row) = wT*standard + wD.*fd;
            jnd_2D2Dra(row) = sqrt( (wT*jndD1/c2D2Dra(subj,1)).^2 + (wD.*jndD2./afs(c2D2Dra(subj,1),c2D2Dra(subj,2),c2D2Dra(subj,3),prop,(standard-fd))).^2 );
            mcp_2D2Dra(row,:) = 0.5 * ( 1+erf( (fc - pse_2D2Dra(row) ) ./ sqrt( 2 * jnd_2D2Dra(row)^2 ) ) );
            
            rsscp_2D2Dra(row) = sum( (mcp_2D2Dra(row,:)-out.choiceprob(row,:)).^2 );
            rsspse_2D2Dra(row) = sum( (pse_2D2Dra(row)-out.pse(row)).^2 );
            rssjnd_2D2Dra(row) = sum( (jnd_2D2Dra(row)-out.jnd(row)).^2 );
            

            wT = c2D1Daa(subj,4);
            wD = afsp(c2D1Daa(subj,4),c2D1Daa(subj,5),prop);
            pse_2D1Daa(row) = wT*standard + wD.*fd;
            jnd_2D1Daa(row) = sqrt( (wT*jndD1/c2D1Daa(subj,1)).^2 + (wD.*jndD2./afs(c2D1Daa(subj,1),c2D1Daa(subj,2),c2D1Daa(subj,3),prop,(standard-fd))).^2 );
            mcp_2D1Daa(row,:) = 0.5 * ( 1+erf( (fc - pse_2D1Daa(row) ) ./ sqrt( 2 * jnd_2D1Daa(row)^2 ) ) );
            
            rsscp_2D1Daa(row) = sum( (mcp_2D1Daa(row,:)-out.choiceprob(row,:)).^2 );
            rsspse_2D1Daa(row) = sum( (pse_2D1Daa(row)-out.pse(row)).^2 );
            rssjnd_2D1Daa(row) = sum( (jnd_2D1Daa(row)-out.jnd(row)).^2 );
            
            
            wT = 1 ./ (c2D1Dra(subj,4) + afspn(c2D1Dra(subj,5),prop));
            wD = afspn(c2D1Dra(subj,5),prop) ./ (c2D1Dra(subj,4) + afspn(c2D1Dra(subj,5),prop));
            pse_2D1Dra(row) = wT*standard + wD.*fd;
            jnd_2D1Dra(row) = sqrt( (wT*jndD1/c2D1Dra(subj,1)).^2 + (wD.*jndD2./afs(c2D1Dra(subj,1),c2D1Dra(subj,2),c2D1Dra(subj,3),prop,(standard-fd))).^2 );
            mcp_2D1Dra(row,:) = 0.5 * ( 1+erf( (fc - pse_2D1Dra(row) ) ./ sqrt( 2 * jnd_2D1Dra(row)^2 ) ) );
            
            rsscp_2D1Dra(row) = sum( (mcp_2D1Dra(row,:)-out.choiceprob(row,:)).^2 );
            rsspse_2D1Dra(row) = sum( (pse_2D1Dra(row)-out.pse(row)).^2 );
            rssjnd_2D1Dra(row) = sum( (jnd_2D1Dra(row)-out.jnd(row)).^2 );
            
            
            wT = c2D1Dpaa(subj,4);
            wD = afsp(c2D1Dpaa(subj,4),c2D1Dpaa(subj,3),prop);
            pse_2D1Dpaa(row) = wT*standard + wD.*fd;
            jnd_2D1Dpaa(row) = sqrt( (wT*jndD1/c2D1Dpaa(subj,1)).^2 + (wD.*jndD2./afs(c2D1Dpaa(subj,1),c2D1Dpaa(subj,2),c2D1Dpaa(subj,3),prop,(standard-fd))).^2 );
            mcp_2D1Dpaa(row,:) = 0.5 * ( 1+erf( (fc - pse_2D1Dpaa(row) ) ./ sqrt( 2 * jnd_2D1Dpaa(row)^2 ) ) );
            
            rsscp_2D1Dpaa(row) = sum( (mcp_2D1Dpaa(row,:)-out.choiceprob(row,:)).^2 );
            rsspse_2D1Dpaa(row) = sum( (pse_2D1Dpaa(row)-out.pse(row)).^2 );
            rssjnd_2D1Dpaa(row) = sum( (jnd_2D1Dpaa(row)-out.jnd(row)).^2 );
            
            
            wT = 1 ./ (c2D1Dpra(subj,4) + afspn(c2D1Dpra(subj,3),prop));
            wD = afspn(c2D1Dpra(subj,3),prop) ./ (c2D1Dpra(subj,4) + afspn(c2D1Dpra(subj,3),prop));
            pse_2D1Dpra(row) = wT*standard + wD.*fd;
            jnd_2D1Dpra(row) = sqrt( (wT*jndD1/c2D1Dpra(subj,1)).^2 + (wD.*jndD2./afs(c2D1Dpra(subj,1),c2D1Dpra(subj,2),c2D1Dpra(subj,3),prop,(standard-fd))).^2 );
            mcp_2D1Dpra(row,:) = 0.5 * ( 1+erf( (fc - pse_2D1Dpra(row) ) ./ sqrt( 2 * jnd_2D1Dpra(row)^2 ) ) );
            
            rsscp_2D1Dpra(row) = sum( (mcp_2D1Dpra(row,:)-out.choiceprob(row,:)).^2 );
            rsspse_2D1Dpra(row) = sum( (pse_2D1Dpra(row)-out.pse(row)).^2 );
            rssjnd_2D1Dpra(row) = sum( (jnd_2D1Dpra(row)-out.jnd(row)).^2 );
            

            
            wT = (jndD2./afs(c2Drr(subj,1),c2Drr(subj,2),c2Drr(subj,3),prop,(standard-fd))).^2 ./ ( (jndD1/c2Drr(subj,1))^2 + (jndD2./afs(c2Drr(subj,1),c2Drr(subj,2),c2Drr(subj,3),prop,(standard-fd))).^2 );
            wD = (jndD1/c2Drr(subj,1))^2 ./ ( (jndD1/c2Drr(subj,1))^2 + (jndD2./afs(c2Drr(subj,1),c2Drr(subj,2),c2Drr(subj,3),prop,(standard-fd))).^2 );
            pse_2Drr(row) = wT*standard + wD.*fd;
            jnd_2Drr(row) = sqrt( wT.^2 .* (jndD2./afs(c2Drr(subj,1),c2Drr(subj,2),c2Drr(subj,3),prop,(standard-fd))).^2 + wD.^2 * (jndD1/c2Drr(subj,1))^2 );
            mcp_2Drr(row,:) = 0.5 * ( 1+erf( (fc - pse_2Drr(row) ) ./ sqrt( 2 * jnd_2Drr(row)^2 ) ) );
            
            rsscp_2Drr(row) = sum( (mcp_2Drr(row,:)-out.choiceprob(row,:)).^2 );
            rsspse_2Drr(row) = sum( (pse_2Drr(row)-out.pse(row)).^2 );
            rssjnd_2Drr(row) = sum( (jnd_2Drr(row)-out.jnd(row)).^2 );
            
            
            wT = c2Daa(subj,1);
            wD = afs(c2Daa(subj,1),c2Daa(subj,2),c2Daa(subj,3),prop,(standard-fd));
            pse_2Daa(row) = wT*standard + wD.*fd;
            jnd_2Daa(row) = sqrt( jndD1^2 + jndD2^2 );
            mcp_2Daa(row,:) = 0.5 * ( 1+erf( (fc - pse_2Daa(row) ) ./ sqrt( 2 * jnd_2Daa(row)^2 ) ) );
            
            rsscp_2Daa(row) = sum( (mcp_2Daa(row,:)-out.choiceprob(row,:)).^2 );
            rsspse_2Daa(row) = sum( (pse_2Daa(row)-out.pse(row)).^2 );
            rssjnd_2Daa(row) = sum( (jnd_2Daa(row)-out.jnd(row)).^2 );
            
            
            wT = 1./(c2Dra(subj,4)+afsn(c2Dra(subj,2),c2Dra(subj,3),prop,(standard-fd)));
            wD = afsn(c2Dra(subj,2),c2Dra(subj,3),prop,(standard-fd))./(c2Dra(subj,4)+afsn(c2Dra(subj,2),c2Dra(subj,3),prop,(standard-fd)));
            pse_2Dra(row) = wT*standard + wD.*fd;
            jnd_2Dra(row) = sqrt( ( wT * jndD1 / c2Dra(subj,1) ).^2 + ( wD * jndD2./afs(c2Dra(subj,1),c2Dra(subj,2),c2Dra(subj,3),prop,(standard-fd)) ).^2 );
            mcp_2Dra(row,:) = 0.5 * ( 1+erf( (fc - pse_2Dra(row) ) ./ sqrt( 2 * jnd_2Dra(row)^2 ) ) );
            
            rsscp_2Dra(row) = sum( (mcp_2Dra(row,:)-out.choiceprob(row,:)).^2 );
            rsspse_2Dra(row) = sum( (pse_2Dra(row)-out.pse(row)).^2 );
            rssjnd_2Dra(row) = sum( (jnd_2Dra(row)-out.jnd(row)).^2 );
            
            
            wT = (jndD2./afsp(c1Dprr(subj,1),c1Dprr(subj,2),prop)).^2 ./ ( (jndD1/c1Dprr(subj,1))^2 + (jndD2./afsp(c1Dprr(subj,1),c1Dprr(subj,2),prop)).^2);
            wD = (jndD1/c1Dprr(subj,1))^2 ./ ( (jndD1/c1Dprr(subj,1))^2 + (jndD2./afsp(c1Dprr(subj,1),c1Dprr(subj,2),prop)).^2);
            pse_1Dprr(row) = wT*standard + wD.*fd;
            jnd_1Dprr(row) = sqrt( (wT.*jndD1/c1Dprr(subj,1)).^2 + (wD.*jndD2./afsp(c1Dprr(subj,1),c1Dprr(subj,2),prop)).^2 );
            mcp_1Dprr(row,:) = 0.5 * ( 1+erf( (fc - pse_1Dprr(row) ) ./ sqrt( 2 * jnd_1Dprr(row)^2 ) ) );
            
            rsscp_1Dprr(row) = sum( (mcp_1Dprr(row,:)-out.choiceprob(row,:)).^2 );
            rsspse_1Dprr(row) = sum( (pse_1Dprr(row)-out.pse(row)).^2 );
            rssjnd_1Dprr(row) = sum( (jnd_1Dprr(row)-out.jnd(row)).^2 );
            
            
            wT = c1Dpaa(1);
            wD = afsp(c1Dpaa(subj,1),c1Dpaa(subj,2),prop);
            pse_1Dpaa(row) = wT*standard + wD.*fd;
            jnd_1Dpaa(row) = sqrt( jndD1^2 + jndD2^2 );
            mcp_1Dpaa(row,:) = 0.5 * ( 1+erf( (fc - pse_1Dpaa(row) ) ./ sqrt( 2 * jnd_1Dpaa(row)^2 ) ) );
            
            rsscp_1Dpaa(row) = sum( (mcp_1Dpaa(row,:)-out.choiceprob(row,:)).^2 );
            rsspse_1Dpaa(row) = sum( (pse_1Dpaa(row)-out.pse(row)).^2 );
            rssjnd_1Dpaa(row) = sum( (jnd_1Dpaa(row)-out.jnd(row)).^2 );
            
            
            wT = c1Dpra(subj,1)./(c1Dpra(subj,3)+c1Dpra(subj,1)+afsp(c1Dpra(subj,1),c1Dpra(subj,2),prop));
            wD = afsp(c1Dpra(subj,1),c1Dpra(subj,2),prop)./(c1Dpra(subj,3)+c1Dpra(subj,1)+afsp(c1Dpra(subj,1),c1Dpra(subj,2),prop));
            pse_1Dpra(row) = wT*standard + wD.*fd;
            jnd_1Dpra(row) = sqrt( (wT.*jndD1/c1Dpra(subj,1)).^2 + (wD.*jndD2./afsp(c1Dpra(subj,1),c1Dpra(subj,2),prop)).^2);
            mcp_1Dpra(row,:) = 0.5 * ( 1+erf( (fc - pse_1Dpra(row) ) ./ sqrt( 2 * jnd_1Dpra(row)^2 ) ) );
            
            rsscp_1Dpra(row) = sum( (mcp_1Dpra(row,:)-out.choiceprob(row,:)).^2 );
            rsspse_1Dpra(row) = sum( (pse_1Dpra(row)-out.pse(row)).^2 );
            rssjnd_1Dpra(row) = sum( (jnd_1Dpra(row)-out.jnd(row)).^2 );
            
            
            wT = (jndD2./afsf(c1Dfrr(subj,1),c1Dfrr(subj,2),(standard-fd))).^2 ./ ( (jndD1/c1Dfrr(subj,1))^2 + (jndD2./afsf(c1Dfrr(subj,1),c1Dfrr(subj,2),(standard-fd))).^2);
            wD = (jndD1/c1Dfrr(subj,1))^2 ./ ( (jndD1/c1Dfrr(subj,1))^2 + (jndD2./afsf(c1Dfrr(subj,1),c1Dfrr(subj,2),(standard-fd))).^2);
            pse_1Dfrr(row) = wT*standard + wD.*fd;
            jnd_1Dfrr(row) = sqrt( (wT.*jndD1/c1Dfrr(subj,1)).^2 + (wD.*jndD2./afsp(c1Dfrr(subj,1),c1Dfrr(subj,2),(standard-fd))).^2);
            mcp_1Dfrr(row,:) = 0.5 * ( 1+erf( (fc - pse_1Dfrr(row) ) ./ sqrt( 2 * jnd_1Dfrr(row)^2 ) ) );
            
            rsscp_1Dfrr(row) = sum( (mcp_1Dfrr(row,:)-out.choiceprob(row,:)).^2 );
            rsspse_1Dfrr(row) = sum( (pse_1Dfrr(row)-out.pse(row)).^2 );
            rssjnd_1Dfrr(row) = sum( (jnd_1Dfrr(row)-out.jnd(row)).^2 );
            
            
            wT = c1Dfaa(1);
            wD = afsf(c1Dfaa(subj,1),c1Dfaa(subj,2),(standard-fd));
            pse_1Dfaa(row) = wT*standard + wD.*fd;
            jnd_1Dfaa(row) = sqrt( jndD1^2 + jndD2^2 );
            mcp_1Dfaa(row,:) = 0.5 * ( 1+erf( (fc - pse_1Dfaa(row) ) ./ sqrt( 2 * jnd_1Dfaa(row)^2 ) ) );
            
            rsscp_1Dfaa(row) = sum( (mcp_1Dfaa(row,:)-out.choiceprob(row,:)).^2 );
            rsspse_1Dfaa(row) = sum( (pse_1Dfaa(row)-out.pse(row)).^2 );
            rssjnd_1Dfaa(row) = sum( (jnd_1Dfaa(row)-out.jnd(row)).^2 );
            
            
            wT = c1Dfra(subj,1)./(c1Dfra(subj,3)+c1Dfra(subj,1)+afsf(c1Dfra(subj,1),c1Dfra(subj,2),(standard-fd)));
            wD = afsf(c1Dfra(subj,1),c1Dfra(subj,2),(standard-fd))./(c1Dfra(subj,3)+c1Dfra(subj,1)+afsf(c1Dfra(subj,1),c1Dfra(subj,2),(standard-fd)));
            pse_1Dfra(row) = wT*standard + wD.*fd;
            jnd_1Dfra(row) = sqrt( (wT.*jndD1/c1Dfra(subj,1)).^2 + (wD.*jndD2./afsf(c1Dfra(subj,1),c1Dfra(subj,2),(standard-fd))).^2);
            mcp_1Dfra(row,:) = 0.5 * ( 1+erf( (fc - pse_1Dfra(row) ) ./ sqrt( 2 * jnd_1Dfra(row)^2 ) ) );
            
            rsscp_1Dfra(row) = sum( (mcp_1Dfra(row,:)-out.choiceprob(row,:)).^2 );
            rsspse_1Dfra(row) = sum( (pse_1Dfra(row)-out.pse(row)).^2 );
            rssjnd_1Dfra(row) = sum( (jnd_1Dfra(row)-out.jnd(row)).^2 );
            
            
            wT = c1Df1Dpaa(subj,3);
            wD = afsp(c1Df1Dpaa(subj,3),c1Df1Dpaa(subj,4),prop);
            pse_DfDpaa(row) = wT*standard + wD.*fd;
            jnd_DfDpaa(row) = sqrt( (wT*jndD1/c1Df1Dpaa(subj,1)).^2 + (wD*jndD2./afsf(c1Df1Dpaa(subj,1),c1Df1Dpaa(subj,2),(standard-fd))).^2 );
            mcp_DfDpaa(row,:) = 0.5 * ( 1+erf( (fc - pse_DfDpaa(row) ) ./ sqrt( 2 * jnd_DfDpaa(row)^2 ) ) );
            
            rsscp_DfDpaa(row) = sum( (mcp_DfDpaa(row,:)-out.choiceprob(row,:)).^2 );
            rsspse_DfDpaa(row) = sum( (pse_DfDpaa(row)-out.pse(row)).^2 );
            rssjnd_DfDpaa(row) = sum( (jnd_DfDpaa(row)-out.jnd(row)).^2 );
            
            
            wT = 1 ./ ( c1Df1Dpra(subj,4)+afspn(c1Df1Dpra(subj,3),prop) );
            wD = afspn(c1Df1Dpra(subj,3),prop) ./ ( c1Df1Dpra(subj,4)+afspn(c1Df1Dpra(subj,3),prop) );
            pse_DfDpra(row) = wT*standard + wD.*fd;
            jnd_DfDpra(row) = sqrt( (wT*jndD1/c1Df1Dpra(subj,1)).^2 + (wD*jndD2./afsf(c1Df1Dpra(subj,1),c1Df1Dpra(subj,2),(standard-fd))).^2 );
            mcp_DfDpra(row,:) = 0.5 * ( 1+erf( (fc - pse_DfDpra(row) ) ./ sqrt( 2 * jnd_DfDpra(row)^2 ) ) );
            
            rsscp_DfDpra(row) = sum( (mcp_DfDpra(row,:)-out.choiceprob(row,:)).^2 );
            rsspse_DfDpra(row) = sum( (pse_DfDpra(row)-out.pse(row)).^2 );
            rssjnd_DfDpra(row) = sum( (jnd_DfDpra(row)-out.jnd(row)).^2 );
            
            
            wT = c1Dp1Dfaa(3);
            wD = afsf(c1Dp1Dfaa(subj,3),c1Dp1Dfaa(subj,4),(standard-fd));
            pse_DpDfaa(row) = wT*standard + wD.*fd;
            jnd_DpDfaa(row) = sqrt( (wT*jndD1/c1Dp1Dfaa(subj,1)).^2 + (wD*jndD2./afsp(c1Dp1Dfaa(subj,1),c1Dp1Dfaa(subj,2),prop)).^2 );
            mcp_DpDfaa(row,:) = 0.5 * ( 1+erf( (fc - pse_DpDfaa(row) ) ./ sqrt( 2 * jnd_DpDfaa(row)^2 ) ) );
            
            rsscp_DpDfaa(row) = sum( (mcp_DpDfaa(row,:)-out.choiceprob(row,:)).^2 );
            rsspse_DpDfaa(row) = sum( (pse_DpDfaa(row)-out.pse(row)).^2 );
            rssjnd_DpDfaa(row) = sum( (jnd_DpDfaa(row)-out.jnd(row)).^2 );
            
            
            wT = 1 ./ ( c1Dp1Dfra(subj,4)+afsfn(c1Dp1Dfra(subj,3),(standard-fd)) );
            wD = afsfn(c1Dp1Dfra(subj,3),(standard-fd)) ./ ( c1Dp1Dfra(subj,4)+afsfn(c1Dp1Dfra(subj,3),(standard-fd)) );
            pse_DpDfra(row) = wT*standard + wD.*fd;
            jnd_DpDfra(row) = sqrt( (wT*jndD1/c1Dp1Dfra(subj,1)).^2 + (wD*jndD2./afsp(c1Dp1Dfra(subj,1),c1Dp1Dfra(subj,2),prop)).^2 );
            mcp_DpDfra(row,:) = 0.5 * ( 1+erf( (fc - pse_DpDfra(row) ) ./ sqrt( 2 * jnd_DpDfra(row)^2 ) ) );
            
            rsscp_DpDfra(row) = sum( (mcp_DpDfra(row,:)-out.choiceprob(row,:)).^2 );
            rsspse_DpDfra(row) = sum( (pse_DpDfra(row)-out.pse(row)).^2 );
            rssjnd_DpDfra(row) = sum( (jnd_DpDfra(row)-out.jnd(row)).^2 );
            
            
            mcp_null(row,:) = 0.5 * ( 1 + erf(...
                (fc - standard) ./ sqrt(2 * jndD1^2)...
                ));
            
            rsscp_null(row) = sum( (mcp_null(row,:)-out.choiceprob(row,:)).^2 );
            rsspse_null(row) = sum( (standard-out.pse(row)).^2 );
            rssjnd_null(row) = sum( (jndD1-out.jnd(row)).^2 );
            
            mcp_free1(row,:) = 0.5 * ( 1+erf( (fc - out.pse(row) ) ./ sqrt( 2 * out.jnd(row)^2 ) ) );
            mcp_free(row,:) = 0.5 * ( 1+erf( (fc - out.pse(row) ) ./ sqrt( 2 * out.jnd(row)^2 ) ) );
            
            errorcp_free(row) = sum( (mcp_free1(row,:)-out.choiceprob(row,:)).^2 );
            errorpse_free(row) = sum( (out.pse(row)-out.pse(row)).^2 );
            errorjnd_free(row) = sum( (out.jnd(row)-out.jnd(row)).^2 );
        end
    end
end

performance.mcp = mcp_DfDpra;
performance.ocp = out.choiceprob;

performance.llfree = llfree;
performance.ll = ll;
performance.aic = aic;
performance.bic = bic;
performance.lrt = lrt;


%
clear rsscp rsspse rssjnd
% rss = sqrt( sum (y - yhat)^2 ) % error in Hz/condition
for subj = 1:8
    indx = 9*subj-8:9*subj;
    % RSS
    rsscp(subj,:) = sqrt( mean( [rsscp_2D2Daa(indx)' rsscp_2D2Dra(indx)' rsscp_2D1Daa(indx)' rsscp_2D1Dra(indx)' rsscp_2D1Dpaa(indx)' rsscp_2D1Dpra(indx)' rsscp_2Drr(indx)' rsscp_2Daa(indx)' rsscp_2Dra(indx)' rsscp_1Dprr(indx)' rsscp_1Dpaa(indx)' rsscp_1Dpra(indx)' rsscp_1Dfrr(indx)' rsscp_1Dfaa(indx)' rsscp_1Dfra(indx)' rsscp_DfDpaa(indx)' rsscp_DfDpra(indx)' rsscp_DpDfaa(indx)' rsscp_DpDfra(indx)'] ));
    rsspse(subj,:) = sqrt( mean( [rsspse_2D2Daa(indx)' rsspse_2D2Dra(indx)' rsspse_2D1Daa(indx)' rsspse_2D1Dra(indx)' rsspse_2D1Dpaa(indx)' rsspse_2D1Dpra(indx)' rsspse_2Drr(indx)' rsspse_2Daa(indx)' rsspse_2Dra(indx)' rsspse_1Dprr(indx)' rsspse_1Dpaa(indx)' rsspse_1Dpra(indx)' rsspse_1Dfrr(indx)' rsspse_1Dfaa(indx)' rsspse_1Dfra(indx)' rsspse_DfDpaa(indx)' rsspse_DfDpra(indx)' rsspse_DpDfaa(indx)' rsspse_DpDfra(indx)'] ));
    rssjnd(subj,:) = sqrt( mean( [rssjnd_2D2Daa(indx)' rssjnd_2D2Dra(indx)' rssjnd_2D1Daa(indx)' rssjnd_2D1Dra(indx)' rssjnd_2D1Dpaa(indx)' rssjnd_2D1Dpra(indx)' rssjnd_2Drr(indx)' rssjnd_2Daa(indx)' rssjnd_2Dra(indx)' rssjnd_1Dprr(indx)' rssjnd_1Dpaa(indx)' rssjnd_1Dpra(indx)' rssjnd_1Dfrr(indx)' rssjnd_1Dfaa(indx)' rssjnd_1Dfra(indx)' rssjnd_DfDpaa(indx)' rssjnd_DfDpra(indx)' rssjnd_DpDfaa(indx)' rssjnd_DpDfra(indx)'] ));
    % RSS null model
    rsscpnull(subj) = mean(rsscp_null(indx));
    rsspsenull(subj) = mean(rsspse_null(indx));
    rssjndnull(subj) = mean(rssjnd_null(indx));
    % RSS noise ceiling
    rsscpns(subj) = sqrt(mean(errorcp_free(indx)));
    rsspsens(subj) = sqrt(mean(errorpse_free(indx)));
    rssjndns(subj) = sqrt(mean(errorjnd_free(indx)));

end

rsstotal = sqrt(rsspse.^2 + rssjnd.^2);

for i = 1:19
    [hh,pp,~,stats] = ttest2(rsscp(:,i), rsscpns);
    p(i) = pp;
    t(i) = stats.tstat;
end
performance.rsscp = rsscp;
performance.rsscp_t = t;
performance.rsscp_p = p;

for i = 1:19
    [hh,pp,~,stats] = ttest2(rsspse(:,i), rsspsens);
    p(i) = pp;
    t(i) = stats.tstat;
end
performance.rsspse = rsspse;
performance.rsspse_t = t;
performance.rsspse_p = p;

for i = 1:19
    [hh,pp,~,stats] = ttest2(rssjnd(:,i), rssjndns);
    p(i) = pp;
    t(i) = stats.tstat;
end
performance.rssjnd = rssjnd;
performance.rssjnd_t = t;
performance.rssjnd_p = p;

% 3 classes of models
n =   [17 2 4 6 9 12 15 19]; % Normalization
non = [16 1 3 5 8 11 14 18]; % No Normalization
oci = [7 10 13]; % Optimal cue integration
model_indx = [n non oci];

% Printing model selection criteria; generates supplementary table 2
clc;
fprintf('LL\tLRT\tAIC\tBIC\tErrCP\tErrPSE\tErrJND\tErrTotal\n')
fprintf(['%0.1f' char(177) '%0.1f\t'...
    '%0.2f' char(177) '%0.2f\t'...
    '%0.1f' char(177) '%0.1f\t'...
    '%0.1f' char(177) '%0.1f\t'...
    '%0.2f' char(177) '%0.2f\t'...
    '%0.1f' char(177) '%0.1f\t'...
    '%0.1f' char(177) '%0.1f\t'...
    '%0.1f' char(177) '%0.1f\t'...
    '\n'],...
    [nanmean(ll(:,model_indx));nanstd(ll(:,model_indx))/sqrt(8);...
    nanmean(lrt(:,model_indx));nanstd(lrt(:,model_indx))/sqrt(8);...
    nanmean(aic(:,model_indx));nanstd(aic(:,model_indx))/sqrt(8);...
    nanmean(bic(:,model_indx));nanstd(bic(:,model_indx))/sqrt(8);...
    nanmean(rsscp(:,model_indx));nanstd(rsscp(:,model_indx))/sqrt(8);...
    nanmean(rsspse(:,model_indx));nanstd(rsspse(:,model_indx))/sqrt(8);...
    nanmean(rssjnd(:,model_indx));nanstd(rssjnd(:,model_indx))/sqrt(8);...
    nanmean(rsstotal(:,model_indx));nanstd(rsstotal(:,model_indx))/sqrt(8);...
    ])



% Plot Akaike weight; generate supplementary figure 3
aw = exp(-(nanmean(aic)' - min(nanmean(aic)) )/2) / sum(exp(-(nanmean(aic)' - min(nanmean(aic)) )/2) );
d = exp(-(bsxfun(@minus, aic', min(aic')) / 2 ));
s = nansum( exp(-(bsxfun(@minus, aic', min(aic')) / 2 )) );

meanaw = nanmean(bsxfun(@rdivide, d,s),2);
semaw = nanstd(bsxfun(@rdivide, d,s),[],2) / sqrt(8);

meanaw = meanaw(model_indx);
semaw = semaw(model_indx);

c = [repmat([1 0 0],8,1); repmat([0 0 1],8,1); repmat([0 1 0],3,1)];
subplot(211)
X = smartbar(meanaw,0.9,1,'',c);
eb = errorbar(X,meanaw(:),semaw(:),'ok','MarkerSize',1e-10); errorbar_tick(eb,0)


ylabel('Probability')
ylim([0 0.60])
set(gca,'YTick',[0 0.60])
box off

text(7,0.55,'Cue comb. with normalization','color','r','FontSize',16)
text(7,0.45,'Cue comb. without normalization','color','b','FontSize',16)
text(7,0.35,'Optimal cue integration','color','g','FontSize',16)
