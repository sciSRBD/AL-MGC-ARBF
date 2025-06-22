function [ss,g,pf_true,pf_RBF,cov_RBF,timei]=ALR_MCS(ss,s)
tic
dim=100;% dimension of problem
ss=lhsdesign(dim,dim);%Size of initial training points
ss=-4+8*ss;
g=true_objfun(ss);
%% MGC
opts = struct('metric','euclidean','center','unbiased','fast',1,'max',dim); % default parameters
[~, ~, ~,~,dcor]=MGC(ss,g,opts);
DIC=dcor;
plot(dcor,'b')
%% Stratification strategy for weight assignment
meanDIC=mean(DIC);
stdDIC=std(DIC);
find(DIC>meanDIC+4*stdDIC)
interval=max(DIC)-min(DIC);
dim_h=2;
layorIndex=cell(dim_h,1);
for j=1:dim_h
    index=find(DIC<=min(DIC)+j*interval/dim_h & DIC>=min(DIC)+(j-1)*interval/dim_h);
    layorIndex{j}=index;
end
%% 
error=[];
addDoE=[];
lob=[1e-4,1e-4,1e-4]';%MGC-ARBF hyper-parameters tune
upb=[100,10 10]';%MGC-ARBF hyper-parameters tune
Pf=[];
s=randn(1e6,dim);
pf_true=sum(true_objfun(s)<0)/size(s,1)
for i=1:500
    i
    RBF_model= RBF_build(ss,g,lob,upb,layorIndex);% build MGC-ARBF model
    if size(s,1)>1e5
        miug=ones(size(s,1),1);sigmag=miug;
        for j=1:size(s,1)/1e5
            [YRbf,miugi,sigmagi] = RBF_predictor(RBF_model,s((j-1)*1e5+1:(j)*1e5,:));% prediction by MGC-ARBF model
            miug((j-1)*1e5+1:(j)*1e5)=miugi;
            sigmag((j-1)*1e5+1:(j)*1e5)=sigmagi;
        end
    else
        [YRbf,miug,sigmag] = RBF_predictor(RBF_model,s);
    end
    sigmag(addDoE)=eps;
    P_fail=normcdf(-(miug)./sigmag); % probability of failure at a point
    P_wsp=normcdf(-abs(miug)./sigmag); % probability of wrong sign prediction at a point
    Pfi=sum(YRbf<0)
    %% CV_g learning function
    E_plus=miug.*normcdf(miug./sigmag)+sigmag.*normpdf(miug./sigmag);
    E_minus=miug.*normcdf(-miug./sigmag)-sigmag.*normpdf(-miug./sigmag);
    gap=abs(E_plus)+abs(E_minus);
    gap=gap./(sigmag.^2+miug.^2-gap.^2).^0.5;    
    [maxEI,addindex]=min(gap)
    addDoE=[addDoE;addindex];
    %% stopping condition
    Nff=sum(miug<0);
    miu_Nf=sum(P_fail);
    var_Nf=sum(P_fail.*(1-P_fail));
    std_Nf=var_Nf.^0.5;
    num_Nf=1e4;
    rand_Nf=randn(num_Nf,1)*std_Nf+miu_Nf;
    
    miu_Nwsp=sum(P_wsp);
    var_Nwsp=sum(P_wsp.*(1-P_wsp));
    std_Nwsp=var_Nwsp.^0.5;
    rand_Nwsp=randn(num_Nf,1)*std_Nwsp+miu_Nwsp;
    error1=abs(rand_Nwsp./rand_Nf);
    error1=sort(error1,'ascend');
    error1=error1(9750)
    N_lamda=5;
    error=[error;error1];
    if size(error,1)>N_lamda
        if sum(error(end-N_lamda+1:end)<0.02)==N_lamda
            break
        end
    end
    %%
    ss=[ss;s(addindex,:)];% add ss into the DoE
    g=[g;true_objfun(s(addindex,:))];%
    Pf=[Pf Pfi];
end
timei=toc
I=(miug<0);
pf_RBF=mean(I);
cov_RBF=std(I)/mean(I)/size(s,1)^0.5;
end