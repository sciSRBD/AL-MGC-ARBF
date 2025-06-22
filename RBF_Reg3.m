function  [f, Dist] = RBF_Reg3(bf_c,theta,x)
%% Thin plate spline
%% x--training points m*d
m=size(x,1);
Sx=x.*x*theta';%o(nd),n*d*d*1=n*1
% St=Sx;%o(md)
Em=ones(m,1);%
Mtheta=theta'*Em';%o(d*m)
Sxm=Sx*Em';
D_Euc=Sxm+Sxm'-2*x*(x'.*Mtheta);%o(ndm)
Dist=D_Euc;
f=(Dist+bf_c^2).*log(sqrt(Dist+bf_c^2));
% % f=(Dist).*log((1+bf_c*Dist));
% f=exp(-Dist*bf_c);
%% Thin plate spline --loop-based method
% m=size(x,1);
% Dist=zeros(m,m);
% for i=1:m
%     Dist(i,:)=(((repmat(x(i,:),m,1)-x).^2)*theta')';
% end
% % Dist = dist(C,S');Dist=(Dist.^2)';
% f=(Dist+bf_c^2).*log((Dist+bf_c^2).^0.5);
end