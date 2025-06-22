function x=utox(u,dist_type,dist_para)
cdf=erf(u/sqrt(2))*0.5+0.5;%标准正态分布转化到源分布
switch dist_type
    case 6%均匀分布
        x=cdf*(dist_para(2)-dist_para(1))+dist_para(1);
    case 1 %正态分布
        x=u*dist_para(2)+dist_para(1);
    case 11%Gumbel极值1分布
        a=sqrt(6)*dist_para(2)/pi;%标准差得到尺寸参数
        b=dist_para(1)-0.5772*a;%得到位置参数
        % a=dist_para(2);b=dist_para(1);
        x=(-log(-log(cdf)))*a+b;
    case 12%log正态分布
        a=log(dist_para(1))-1/2*log(dist_para(2)/dist_para(1)^2+1);%方差得到尺寸参数
        b=log(dist_para(2)/dist_para(1)^2+1); %得到位置参数
        b=b^0.5;
% a=dist_para(2);b=dist_para(1);
        x=exp(u*b+a);
end
end