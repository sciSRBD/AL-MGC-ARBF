function  [YRbf,miug,sigmag,YRbf_LOO] = RBF_predictor(RBF_model,x)
mS=RBF_model.data.mS;sS=RBF_model.data.sS;
num_x=size(x,1);
xx=(x - repmat(mS,num_x,1)) ./ repmat(sS,num_x,1);
XX=RBF_model.data.XX;
YY=RBF_model.data.YY;
bf_c=RBF_model.bf_c;
theta=RBF_model.theta;
[f_x, ~] = RBF_Reg2(XX,xx,bf_c,theta);

% [f, ~] = RBF_Reg2(XX,XX,bf_c,theta);

f=RBF_model.f;inv_f=RBF_model.inv_f;sY=RBF_model.data.sY;mY=RBF_model.data.mY;
A=f;A_inv=inv_f;M=size(A,1);

beta=f\YY;
YRbf=f_x*beta;
YRbf=YRbf*sY+mY;
YRbf_LOO=zeros(size(xx,1),size(XX,1));

beta_add0=zeros(M,M);
for ii=1:size(XX,1)%o(M^3)
    index=[1:ii-1 ii+1:M];
    a=f(index,index);
    YY_LOO=YY(index);
    beta_LOO=a\YY_LOO;%% o(M^2)
    beta_add0(index,ii)=beta_LOO;%
end
YRbf_LOO=f_x*beta_add0;% o(M^2*N)
YRbf_LOO=YRbf_LOO*sY+mY;
%     pseudo_YRbf=n*repmat(YRbf,1,n)-(n-1)*YRbf_LOO;
miug=mean(YRbf_LOO,2);% o(M*N)
sigmag=std(YRbf_LOO,0,2);% o(M*N)
end