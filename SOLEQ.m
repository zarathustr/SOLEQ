function q=SOLEQ(Db,Dr,k)

n=length(Db(1,:));

if n>1

    [V1,S]=eig((W_matrix(Db(:,1),Dr(:,1))+eye(4))*0.5);
    [Vn,S]=eig((W_matrix(Db(:,n),Dr(:,n))+eye(4))*0.5);
    S=diag([0 0 1 1]);

%     V1=V_matrix(Db(:,1),Dr(:,1));
%     Vn=V_matrix(Db(:,n),Dr(:,n));

    if((n-2)>0)
        P=eye(4);
        for i=2:n-1
            P=P*(W_matrix(Db(:,i),Dr(:,i))+eye(4))*0.5;
        end
        N=[
            0,0,0,0;
            0,0,0,0;
            V1(1,3),V1(2,3),V1(3,3),V1(4,3);
            V1(1,4),V1(2,4),V1(3,4),V1(4,4);
        ] ...
        *P* ...
        [
            0,0,Vn(1,3),Vn(1,4);
            0,0,Vn(2,3),Vn(2,4);
            0,0,Vn(3,3),Vn(3,4);
            0,0,Vn(4,3),Vn(4,4);
        ];
        
        n1=N(3,3);
        n2=N(3,4);
        n3=N(4,3);
        n4=N(4,4);
    else
        n1=V1(1,3)*Vn(1,3) + V1(2,3)*Vn(2,3) + V1(3,3)*Vn(3,3) + V1(4,3)*Vn(4,3);
        n2=V1(1,3)*Vn(1,4) + V1(2,3)*Vn(2,4) + V1(3,3)*Vn(3,4) + V1(4,3)*Vn(4,4);
        n3=V1(1,4)*Vn(1,3) + V1(2,4)*Vn(2,3) + V1(3,4)*Vn(3,3) + V1(4,4)*Vn(4,3);
        n4=V1(1,4)*Vn(1,4) + V1(2,4)*Vn(2,4) + V1(3,4)*Vn(3,4) + V1(4,4)*Vn(4,4);
    end

    u1=n1*n1 + n2*n2;
    u2=n1*n3 + n2*n4;
    u4=n3*n3 + n4*n4;

    sqrt_delta=sqrt(u1*u1-2*u1*u4 +4*u2*u2 +u4*u4);

    lambda_u1=(-sqrt_delta+u1+u4)*0.5;
    lambda_u2=(sqrt_delta+u1+u4)*0.5;


    VU=[1,0,0,0;
        0,1,0,0;
        0,0,(-u1+u4+sqrt_delta)/u2*0.5,(u1-u4+sqrt_delta)/u2*0.5;
        0,0,1,1];

    VV=V1*VU;

    q=VV*diag([0,0,lambda_u1^k,lambda_u2^k])*VV'*[1;0;0;0];
    q=q./norm(q);
else
    q=((W_matrix(Db(:,1),Dr(:,1))+eye(4))*0.5)*[1;0;0;0];
    q=q./norm(q);
end

end