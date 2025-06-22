function  [f, Dist] = RBF_Reg2(S,x,bf_c,theta)
%% Thin plate spline
% x--test points mx*d
%S--training points m*d

m=size(S,1);
mx=size(x,1);
Sx=x.*x*theta';%o(nd),n*d*d*1=n*1
St=S.*S*theta';%o(md)
Em=ones(m,1);%size of trianing points
En=ones(mx,1);%size of test points
Mtheta=theta'*Em';%o(d*m)
D_Euc=Sx*Em'+En*St'-2*x*(S'.*Mtheta);%o(ndm)
Dist=D_Euc;
f=(Dist+bf_c^2).*log(sqrt(Dist+bf_c^2));% Thin Plate spline
%% Thin plate spline --loop-based method
% m=size(x,1);
% Dist=zeros(m,m);
% for i=1:m
%     Dist(i,:)=(((repmat(x(i,:),m,1)-x).^2)*theta')';
% end
% % Dist = dist(C,S');Dist=(Dist.^2)';
% f=(Dist+bf_c^2).*log((Dist+bf_c^2).^0.5);
end