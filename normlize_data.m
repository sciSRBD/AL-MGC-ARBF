function [data]=normlize_data(X,Y)
num_ss=size(X,1);
mS = mean(X);   sS = std(X);
mY = mean(Y);   sY = std(Y);
XX = (X - repmat(mS,num_ss,1)) ./ repmat(sS,num_ss,1);
YY = (Y - repmat(mY,num_ss,1)) ./ repmat(sY,num_ss,1);
data.XX=XX;data.YY=YY;data.mS=mS;data.sS=sS;data.mY=mY;data.sY=sY;
end