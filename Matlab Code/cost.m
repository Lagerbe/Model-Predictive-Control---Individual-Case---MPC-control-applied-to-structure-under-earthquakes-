function J=cost(init,X1,u,RHO,RHO2,Nu,Nu2,Nd,Nd2,Dd,Cy,Cz,Du,xg,v,Ne1,Ne2,e,P12,P1,Q,R,N)

    [X1,Y,Z]=NL_system(init,X1,u(1),RHO,RHO2,Nu,Nu2,Nd,Nd2,Dd,Cy,Cz,Du,0,0,Ne1,Ne2,e,P12,P1,1);
    J=0;
    for k=2:N
        J=J+X1'*Q*X1+u(k)'*R*u(k);
        [X1,Y,Z]=NL_system(init,X1,u(k),RHO,RHO2,Nu,Nu2,Nd,Nd2,Dd,Cy,Cz,Du,0,0,Ne1,Ne2,e,P12,P1,k);
    end