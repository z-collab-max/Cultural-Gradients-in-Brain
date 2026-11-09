function B=roi_midpoint_segments(vertices,faces,mask,varargin)
% Intersect each boundary triangle with the binary ROI midpoint contour.
mask=logical(mask(:));m=mask(faces);crossing=any(m,2)&~all(m,2);f=faces(crossing,:);m=m(crossing,:);
n=size(f,1);a=nan(n,3);b=a;count=zeros(n,1);
for edge=[1 2;2 3;3 1]'
 i=edge(1);j=edge(2);ix=find(m(:,i)~=m(:,j));
 mid=(vertices(f(ix,i),:)+vertices(f(ix,j),:))/2;
 first=count(ix)==0;a(ix(first),:)=mid(first,:);b(ix(~first),:)=mid(~first,:);count(ix)=count(ix)+1;
end
assert(all(count==2));v=nan(n*3,3);v(1:3:end,:)=a;v(2:3:end,:)=b;B={v};
end
