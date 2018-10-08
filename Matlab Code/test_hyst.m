%test of our hysteresis function 
close all
clear all
t=0:0.05:30;
x=0.0018*sin(20*t);
hyst=[];
K1=50;
force=[];
for j=1:size(x,2)
    [K,y1,y2]=hysteresis([],x(j),100*10^9,10,64*10^6,K1);
    hyst=[hyst,y1];
    force=[force,y2];
end
subplot(3,1,1)
plot(x,hyst)
% axis([-0.004 0.004 40000 210000])
title("value of the transition force between elasticity and plasticity")
subplot(3,1,2)
plot(x,force)
title("value of the force respond of the system")
subplot(3,1,3)
plot(x)
title("displacement")