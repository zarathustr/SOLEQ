function [v] = GS4(A)
v(:,1)=A(:,1)/norm(A(:,1));
Num=0;
for j=2:4
    sum(:,1)=[0,0,0,0];
    for i=1:j-1
        hij=A(:,j)'*v(:,i);
        sum(:,1) =sum(:,1)+hij*v(:,i);
        Num=Num+1;
    end
    v(:,j)=A(:,j)-sum(:,1);
    v(:,j)=v(:,j)/sqrt(v(1,j)*v(1,j)+v(2,j)*v(2,j)+v(3,j)*v(3,j)+v(4,j)*v(4,j));
end

end