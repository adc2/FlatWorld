function dd=sqdist(p1,p2)
%dd=sqdist(p1,p2) - squared distance between two arrays of 2D points
%
% p1: array of size n1 X 2
% p2: array of size n2 X 2
%

a=size(p1,1);
b=size(p2,1);
p1=repmat(p1,b,1);
p2=repmat(p2,a,1); 
d1=(reshape(p1(:,1),a,b)-reshape(p2(:,1),b,a)').^2;
d2=(reshape(p1(:,2),a,b)-reshape(p2(:,2),b,a)').^2;
dd=d1+d2;
end % of function sqdist
