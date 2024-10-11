function [W,F] = RBF(points0, f0, points, s, Method)
% 
[m1,n1]=size(points0);
[m2,n2]=size(f0);
[m3,n3]=size(points);
if m1~=m2
    warning("采样点与属性数目不一致！")
end
if n1~=n3
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
DisMat=squareform(pdist(points0));
K3=fun(DisMat, s);
for i=1:n2
    W(:,i)=K3\f0(:,i);
end


%插值
F=zeros(m3,n2);
R=pdist2(points,points0);
R_K=fun(R,s);
f=R_K*W;
F=F+f;

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
