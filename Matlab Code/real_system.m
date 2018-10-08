function [X,Y,Z]=real_system(X1,u,RHO,Nu,Nd,Cz,Dd,Cy,Du,xg,v)
    
    X=RHO*X1+Nu*u+Nd*xg;
    Z=Cz*X1;
    Y=Cy*X1+Du*u+Dd*xg+v;
    
    
    
    