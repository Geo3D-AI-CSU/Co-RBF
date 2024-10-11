function F = idw(points0, f0, points, k)
%UNTITLED2 此处显示有关此函数的摘要
%   此处显示详细说明
D=pdist2(points0,points);
w=D.^(-k);
W=w./sum(w,1);
W(isnan(W))=1;
F=sum(W.*f0,1)';

end