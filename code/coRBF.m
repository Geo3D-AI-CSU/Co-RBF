function [W,F] = coRBF(MVpoints, MVf, CVpoints, CVf, points, r, Method)
% MV为主变量
% CV为主变量
% points为待插点坐标
% Method为基函数方法
% r为基函数作用半径，仅高斯基函数需要使用
% 输出为n+2个参数值

[m1,n1]=size(MVpoints);
[m2,n2]=size(MVf);
[m3,n3]=size(points);
[m4,n4]=size(CVpoints);
[m5,n5]=size(CVf);

% k=10;
% points(:,3)=points(:,3)*k;
% MVpoints(:,3)=MVpoints(:,3)*k; 
% CVpoints(:,3)=CVpoints(:,3)*k;

if m1~=m2
    warning("主变量采样点与属性数目不一致！")
end
if m4~=m5
    warning("次变量采样点与属性数目不一致！")
end
if n1~=n3||n1~=n4||n4~=n3
    warning("待插点与采样点维度不一致！")
end

%选择基函数
switch Method
    case 'linear'
        fun=@Kernel_Linear;
    case 'gaussian'
        fun=@Kernel_Gaussian;
    case 'cubic'
        fun=@Kernel_Cubic;
    case 'thin_plate'
        fun=@Kernel_Thin_plate;
end


%计算权重系数
D1 = pdist2(MVpoints, MVpoints);
A1 = fun(D1,r);
[~,P] = RBF(CVpoints, CVf, MVpoints, r, Method);
A=[A1,P, ones(m1,1);
      P',  0,    0;
   ones(1,m1),0,0];
b=[MVf;0;0];
W=A\b;
% A=[A1,P;
%       P',  0;
%    ones(1,m1),0];
% b=[MVf;0;0];
% W=A\b;


% 插值
D1_K=pdist2(points, MVpoints);
A1_K=fun(D1_K, r);
[~,P] = RBF(CVpoints, CVf, points, r, Method);
F=[A1_K,P,ones(size(points,1),1)]*W;
% F=[A1_K,P]*W;

end


function z=Kernel_Linear(R, s)
    z=abs(R);
end

function z=Kernel_Gaussian(R, s)
    z=exp(-R.^2/2/s^2);
end

function z=Kernel_Cubic(R, s)
    z=R.^3;
end

function z=Kernel_Thin_plate(R, s)
    z=R.^2.*log(R);
end

function P=Potential(points1, points2, f)
    D2=pdist2(points1, points2);
    M=D2.^(-2);
    total_M=sum(M,2);
    w=M./total_M;
    w(isnan(w))=1;
    P=sum(w.*(f'), 2);
end
