function x=utox(u,dist_type,dist_para)
cdf=erf(u/sqrt(2))*0.5+0.5;%��׼��̬�ֲ�ת����Դ�ֲ�
switch dist_type
    case 6%���ȷֲ�
        x=cdf*(dist_para(2)-dist_para(1))+dist_para(1);
    case 1 %��̬�ֲ�
        x=u*dist_para(2)+dist_para(1);
    case 11%Gumbel��ֵ1�ֲ�
        a=sqrt(6)*dist_para(2)/pi;%��׼��õ��ߴ����
        b=dist_para(1)-0.5772*a;%�õ�λ�ò���
        % a=dist_para(2);b=dist_para(1);
        x=(-log(-log(cdf)))*a+b;
    case 12%log��̬�ֲ�
        a=log(dist_para(1))-1/2*log(dist_para(2)/dist_para(1)^2+1);%����õ��ߴ����
        b=log(dist_para(2)/dist_para(1)^2+1); %�õ�λ�ò���
        b=b^0.5;
% a=dist_para(2);b=dist_para(1);
        x=exp(u*b+a);
end
end