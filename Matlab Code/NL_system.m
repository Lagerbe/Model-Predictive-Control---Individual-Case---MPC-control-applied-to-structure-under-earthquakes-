function [X,Y,Z]=NL_system(init,X1,u,RHO,RHO2,Nu,Nu2,Nd,Nd2,Dd,Cy,Cz,Du,xg,v,Ne1,Ne2,e,P12,P1,k)
    global next_pred x_prev x_out_of_deform ed yy x_in
    K1=1387709.943;
    dt=0.02;

    [K,y1,y2]=hysteresis(init,X1(1),90*10^9,K1*96/100,20*10^6,K1);
%     else
%        [K,y1,y2]=hysteresis([],X1(1),35*10^9,K1*98/100,60*10^6,K1);
%     end
%     if k==2
%         init=[ed,x_out_of_deform,x_prev,yy,x_in];
%     end
%         [K,y1,y2]=hysteresis(init,X1(1),35*10^9,K1*97/100,60*10^6,K1);
%     end

    if K>=K1
        Ne=Ne1;
        RHO3=RHO;
        Nu3=Nu;
        Nd3=Nd;
        a=P12*[0;y1*dt];
        b=y1;
        X=RHO3*X1+Nu3*u+Ne*e+Nd3*xg-P1*[0;y1*dt];
        Z=Cz*X1;
        Y=Cy*X1+Du*u+Dd*xg+v;
    else
        Ne=Ne2;
        RHO3=RHO2;
        Nu3=Nu2;
        Nd3=Nd2;
        c=P12*[0;y1*dt];
        d=y1;
        X=RHO3*X1+Nu3*u+Ne*e+Nd3*xg-P12*[0;y1*dt];
        Z=Cz*X1;
        Y=Cy*X1+Du*u+Dd*xg+v;
    end


        