function [Z] = kriging(points0, f0, points, model, h)
% points0,f0为采样点的坐标与属性值
%  points为待插点
% model为拟合变差函数所用模型
% h为滞后距


[m1,n1]=size(points0);
[m2,n2]=size(f0);
[m3,n3]=size(points);
if m1~=m2
    warning("采样点与属性数目不一致！")
end
if n1~=n3
    warning("待插点与采样点维度不一致！")
end

%计算距离矩阵与半方差
D=pdist2(points0,points0);
R=pdist2(f0,f0);
R=0.5*R.^2;
n=ceil(max(max(D))/h);
% dh=ceil(max(D,all)/n);
% 
% for i=1:n
%     r(i)=mean(R(D<i*h&D>=(i-1)*h));
%     d(i)=mean(D(D<i*h&D>=(i-1)*h));
% end
% d_r_divided=[d',r'];

d_r=[D(:),R(:)];
[idx,~]=kmeans(d_r(:,1),n);
d_r_divided=[];
for i=1:n
    d_r_divided=[d_r_divided; mean(d_r(find(idx==i),:),1)];
end

d_r_divided=sortrows(d_r_divided,1);
figure;
scatter(d_r_divided(:,1),d_r_divided(:,2))

%选择变差函数模型
switch model
    % case 'Pure Nugget Effect'
    case 'Exponential'
        ft=fittype('c0+c*(1-exp(-(3*x./a)))');
        fun=@Var_Exponential;
    case 'Spherical'
        ft=fittype('c0+c*(3*x./2/a-x.^3/2/(a^3))');
        fun=@Var_Spherical;
    case 'Gaussian'
        ft=fittype('c0+c*(1-exp(-(3*x.^2./a^2)))');
            fun=@Var_Gaussian;
end

%拟合变差函数
[fitresult, gof] = fit(d_r_divided(:,1), d_r_divided(:,2), ft);
c0=fitresult.c0;
c=fitresult.c;
a=fitresult.a; %变程
y=fun(d_r_divided(:,1),c0,c,a);
hold on 
x=d_r_divided(:,1);
plot(x,y)

%计算权重，插值
A=[fun(D,c0,c,a ...
    ),ones(m1,1);ones(1,m1),0];
D_K=pdist2(points0,points);
R_K=fun(D_K,c0,c,a);
B=[R_K;ones(1,m3)];
W=A\B;
temp=W(1:m1,:).*f0;
Z=sum(temp,1)';

end

function y=Var_Spherical(x,c0,c,a)
    y=c0+c*(3*x./2/a-x.^3/2/(a^3));
    y(x>a)=c0+c;
end

function y=Var_Exponential(x,c0,c,a)
    y=c0+c*(1-exp(-(3*x./a)));
    y(x>a)=c0+c;
end

function y=Var_Gaussian(x,c0,c,a)
    y=c0+c*(1-exp(-(3*x.^2./a^2)));

    y(x>a)=c0+c;
end

