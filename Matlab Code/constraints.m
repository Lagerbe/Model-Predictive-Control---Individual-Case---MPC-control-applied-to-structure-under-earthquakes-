function [Cin,Ceq]=constraints(init,X1,u,RHO,RHO2,Nu,Nu2,Nd,Nd2,Dd,Cy,Cz,Du,xg,v,Ne1,Ne2,e,P12,P1,N)
    Ceq=[];
    X2=[];
    for k=1:N
        [X,Y,Z]=NL_system(init,X1,u,RHO,RHO2,Nu,Nu2,Nd,Nd2,Dd,Cy,Cz,Du,xg,v,Ne1,Ne2,e,P12,P1,k);
        X2=[X2,X];
    end
    Cin(1)=max(X2(1,:))-18/90000;
    Cin(2)=-18/90000-min(X2(1,:));