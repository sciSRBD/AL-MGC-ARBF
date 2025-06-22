function  RBF_model= RBF_build(X,Y,lob,upb,layorIndex)
[data]=normlize_data(X,Y);
XX=data.XX;
YY=data.YY;
nvars=size(layorIndex,1)+1; % dim of hyperrameters
% lob=1e-4*ones(nvars,1);
% upb=100*ones(nvars,1);
%%
popsize=nvars*10; % population
maxgen=100; % iteraion
tolfun1=1e-6;
% Start with the default options
options = gaoptimset;
% Modify options setting
options = gaoptimset(options,'Vectorized','on');
options = gaoptimset(options,'Display','iter');
options = gaoptimset(options,'PopulationType', 'doubleVector');
options = gaoptimset(options,'PopulationSize',popsize);
options = gaoptimset(options,'CreationFcn', @gacreationuniform);
options = gaoptimset(options,'FitnessScalingFcn', @fitscalingrank);
options = gaoptimset(options,'SelectionFcn', @selectionstochunif);
options = gaoptimset(options,'MutationFcn', @mutationadaptfeasible);
% options = gaoptimset(options,'CrossoverFraction', 0.4);
% options = gaoptimset(options,'MigrationFraction', 0.1);
options = gaoptimset(options,'Generations', maxgen);
options = gaoptimset(options,'TolFun',tolfun1);
options.UseParallel='always';
options.Vectorized='on';
% hybridopts = optimset('display','iter','LargeScale','off','Algorithm','SQP');
hybridopts = optimset('display','iter','LargeScale','off','Algorithm','SQP','TolFun',1e-16,'TolX',1e-9);
options = gaoptimset(options,'HybridFcn',{@fmincon,hybridopts});
RBF_parameter=(lob+upb)/2;
f0=shape_tune(RBF_parameter',layorIndex,data);
flag=1;
while flag % multiple tests of optimizaiton
    [RBF_parameter, f]=ga(@(RBF_parameter)shape_tune(RBF_parameter,layorIndex,data),nvars,[],[],[],[],lob',upb',[],options)
    if abs(f0-f)>1e-6
        flag=0
    end
end
%%
theta=RBF_parameter(:,2:end);
bf_c=RBF_parameter(:,1);
thetaA=ones(size(theta,1),size(data.XX,2));
for j=1:size(layorIndex,1)
    Index=layorIndex{j}';
    thetaA(:,Index)=theta(:,j)*ones(1,length(Index));
end
theta=thetaA;
[f, Dist] = RBF_Reg3(bf_c,theta,XX);
if rcond(f)<1e-16
    mu=1e-10;
    f = f + mu*eye(size(f));
    rcond(f);
end
inv_f=f\eye(size(f));
beta=inv_f*YY;
RBF_model.beta=beta;
RBF_model.inv_f=inv_f;
RBF_model.f=f;
RBF_model.data=data;
RBF_model.bf_c=bf_c;
RBF_model.theta=theta;
end
function LOOCV=shape_tune(RBF_parameter,layorIndex,data)
YY=data.YY;
theta=RBF_parameter(:,2:end);% RBF_parameter shape c; layer 1; layer 2
thetaA=ones(size(theta,1),size(data.XX,2));
for j=1:size(layorIndex,1)
    Index=layorIndex{j}';
    thetaA(:,Index)=theta(:,j)*ones(1,length(Index));
end
theta=thetaA;
RBF_c=RBF_parameter(:,1);
LOOCV=zeros(size(RBF_parameter,1),1);
for i=1:size(RBF_parameter,1)
    [f, Dist] = RBF_Reg3(RBF_c(i),theta(i,:),data.XX);    
    if rcond(f)<1e-16
        mu=1e-10;
        f = f + mu*eye(size(f));
        rcond(f);
    end    
    inv_f=f\eye(size(f));
    beta=inv_f*YY;
    error_i=beta./diag(inv_f);
    LOOCV(i)=sum(error_i.^2);
end
end