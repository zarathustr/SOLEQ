function V=V_matrix(Db,Dr)

V=zeros(4,4);
V(1,1)= (Db(1)-Dr(1))*(Db(3)-Dr(3));
V(1,2)=-(Db(1)-Dr(1))*(Db(2)-Dr(2));
V(1,3)= (Db(1)+Dr(1))*(Db(3)+Dr(3));
V(1,4)=-(Db(1)+Dr(1))*(Db(2)+Dr(2));
V(2,1)= (Db(1)-Dr(1))*(Db(2)+Dr(2));
V(2,2)= (Db(1)-Dr(1))*(Db(3)+Dr(3));
V(2,3)= (Db(1)+Dr(1))*(Db(2)-Dr(2));
V(2,4)= (Db(1)+Dr(1))*(Db(3)-Dr(3));
V(3,1)=(Db(2)*Db(2)+Db(3)*Db(3)-Dr(2)*Dr(2)-Dr(3)*Dr(3));
V(3,3)=(Db(2)*Db(2)+Db(3)*Db(3)-Dr(2)*Dr(2)-Dr(3)*Dr(3));
V(4,2)=(Db(2)*Db(2)+Db(3)*Db(3)-Dr(2)*Dr(2)-Dr(3)*Dr(3));
V(4,4)=(Db(2)*Db(2)+Db(3)*Db(3)-Dr(2)*Dr(2)-Dr(3)*Dr(3));


% [Q,R]=qr(V);
% V=Q;

V=GS4(V);

end