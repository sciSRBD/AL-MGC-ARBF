function g=true_objfun(u)
%% A nonlinear function 
x=u;
g=1.5-x(:,end)+0.01*sum(x(:,1:end-1).^2,2);
%% Section 5.2 d=60
% a=830000;
% x=utox(u,1,[3.41,0.2]);
% [N, d] = size(x);
% g = zeros(N,1);  % 初始化输出
% for n = 1:N
%     x_row = x(n, :);  % 第 n 个样本
%     term1 = (x_row(1) - 1)^2;
% 
%     sum_term = 0;
%     for i = 2:d
%         sum_term = sum_term + i * (2 * x_row(i)^2 - x_row(i-1))^2;
%     end
%     g(n) = a - term1 - sum_term;
% end
%% A cantilever beam with 215 random variables
% Sy=utox(u(:,1),1,[600,60])*1e6;%MPa
% w=utox(u(:,2),1,[0.2,0.001]);%m
% h=utox(u(:,3),1,[0.4,0.001]);%m
% Fi=zeros(size(u,1),106);%KN
% lFi=Fi;%m
% for i=1:6
%     Fi(:,i)=utox(u(:,i+3),12,[30+5*i,(2.4+0.4*i)^2]);
% end
% for i=1:6
%     lFi(:,i)=utox(u(:,i+9),1,[4.3+0.1*i,0.01]);
% end
% for i=7:106
%     Fi(:,i)=utox(u(:,i+15-7+1),1,[10,1]);
% end
% for i=7:106
%     lFi(:,i)=utox(u(:,i+115-7+1),1,[0.02*(i-7+1),0.01]);
% end
% Fi=1000*Fi;
% g=Sy-6*sum(Fi.*lFi,2)./w./h.^2;
%% eg Section 5.3 A quadratic function 
% x=u;
% beta=2;
% kappa=5;
% g=beta+(kappa/4)*(x(:,1)-x(:,2)).^2 - sum(x,2)/sqrt(size(x,2));
%% A geodesic space-truss dome
% x=u;
% mean=[2.5e-3*ones(1,5) 2e-3*ones(1,5) 1e-3*ones(1,15) 1.2e-3*ones(1,25) 2.2e-3*ones(1,10) 1.5e-3*ones(1,15) 70 80];
% CV=[0.15*ones(1,5) 0.12*ones(1,5) 0.08*ones(1,15) 0.08*ones(1,25) 0.1*ones(1,10) 0.1*ones(1,15) 0.15 0.15];
% 
% std=mean.*CV;
% for i=1:size(x,2)
%     x(:,i)=utox(u(:,i),1,[mean(i),std(i)]);
% end
% num=size(x,1);
% disp=zeros(size(x,1),1);
% for j=1:num
%     disp(j)=FEM75bar(x(j,:));
% end
% g=33-disp;
end
