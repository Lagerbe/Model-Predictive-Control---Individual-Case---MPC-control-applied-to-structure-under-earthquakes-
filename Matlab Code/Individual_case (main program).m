%% Individual Case
clear all
close all
global next_pred x_prev x_out_of_deform ed yy init x_in
cas=1;

% Fundamental variables
%% Case number of floor equal to 1
if cas==1
    M=2922.7;
    K=[(21.79^2)*M];
    C=[(2*0.0124*21.79)*M];
    G=[zeros(cas,1);-ones(cas,1)];
    kc=371950.8;
    theta=36*pi/180;

    L=-4*kc*cos(theta);

    A=[zeros(1) eye(1);
        -M^(-1)*K -M^(-1)*C];
    B=[zeros(1); M^(-1)*L];
    Cz=eye(2);
    Dz=zeros(2,1);


end
%% Case number of floor equal to 2
if cas==2
    M=974*eye(2);
    K=[2.74 -1.64;
        -1.64 3.02 ]*10^6;
    C=[382.65 -57.27;
        -57.27 456.73];
    G=[zeros(2,1);-ones(2,1)];
    kc=371950.8;
    theta=pi/6;

    L=-4*kc*cos(theta);  

    A=[zeros(2) eye(2);
        -M^(-1)*K -M^(-1)*C];
    B=[zeros(cas); M^(-1)*L];
    Cz=eye(4);
    Dz=zeros(4,1);


end

%% Case number of floor equal to 3 (prepared but not done)

if cas==3
    M=974*eye(3);
    K=[2.74 -1.64 0.37;
        -1.64 3.02 -1.62;
        0.37 -1.62 1.33]*10^6;
    C=[382.65 -57.27 61.64;
        -57.27 456.73 -2.63;
        61.64 -2.63 437.29];
    G=[zeros(3,1);-ones(3,1)];
    kc=372100;
    theta=36*pi/180;

    L=[-4*kc*cos(theta);
        0;
        0];

    A=[zeros(3) eye(3);
        -M^(-1)*K -M^(-1)*C];
    B=[zeros(3); M^(-1)*L];
    Cz=eye(6);
    
    Dz=zeros(6,1);

end
%%
%Sampling time
dt=0.02;

%Horizon 
N=2;

% Matrices of the system

Dd=zeros(cas,1);
RHO=expm(A*dt)

fun=@(t)expm(A.*(t));
P1=integral(fun,0,dt,'ArrayValued',true);

Nu=P1*B;
Nd=P1*G;

Du=-M^(-1)*L;
Cy=[-M^(-1)*K -M^(-1)*C];

% Importation of the model of earthquake
ground=csvread('acceleration_NS.csv',0,1);
ground2=csvread('elcentro_EW.csv',0,1);
ground3=csvread('elcentro_UP.csv',0,1);
t1=(1:size(ground,1))*0.02;
t2=(1:size(ground2,1))*0.02;
t3=(1:size(ground3,1))*0.02;
%% Plot the ground acceleration
% figure
% plot(t1,ground)
% figure
% plot(t2,ground2)
% figure
% plot(t3,ground3)
%%
% Variance of the ground
W=mean((ground-mean(ground)).^2)

% Noise of the measurement  
v=10/769*randn(cas);
V=v*v';



% Initial conditions
x0=zeros(2*cas,1);
u=0;




%% YALMIP PART

% Parameters of the Yalmip and variables
yalmip('clear')

% Number of iteration 
it=600;

%number of state
Nx=2*cas;
%number of input
Nin=1;


% constrain matrices
Q=1*eye(2*cas);
R=200*eye(cas);
%Cm=[eye(cas),zeros(cas)];

%% Comment if we want non linear
% variable of the YALMIP
u=sdpvar(repmat(Nin,1,N),repmat(1,1,N));
x=sdpvar(repmat(Nx,1,N+1),repmat(1,1,N+1));
y=sdpvar(repmat(cas,1,N+1),repmat(1,1,N+1));
%%
X=x0;
X1=x0;
% Initialisation of the arrays for plotting
state(1:2*cas,1)=x0;
e=zeros(cas,1);
error=[];
v2=10/769*randn(cas,1);

P=zeros(2*cas);
for k=1:100
        P=RHO*(P-P*Cy'*((Cy*P*Cy'+V)^(-1))*Cy*P)*RHO'+Nd*W*Nd';
end
Ne1=P*Cy'*((Cy*P*Cy'+V)^(-1));
Ne=P*Cy'*((Cy*P*Cy'+V)^(-1));
%% Plasticity part(only for case 1)
K1=K*96/100;
A=[zeros(cas) eye(cas);
        -M^(-1)*K1 -M^(-1)*C];

    
RHO2=expm(A*dt)

fun=@(t)expm(A.*(t));
P12=integral(fun,0,dt,'ArrayValued',true);

Nu2=P12*B;
Nd2=P12*G;
P=zeros(2*cas)
for k=1:100
        P=RHO*(P-P*Cy'*((Cy*P*Cy'+V)^(-1))*Cy*P)*RHO'+Nd*W*Nd';
end
Ne2=P*Cy'*((Cy*P*Cy'+V)^(-1));
%%
%Computational tima
com_Time=[]

% Preparation of the matrices
energy=[];
input=[];
%% LINEAR CONTROLLER
for j=1:it
    t=cputime;
%     %Noise
%     v=1/769*randn(cas,N);
    Cost_function=0;
    constrain=[x{1}==x0];
    
%     v_cumu=v(:,1);
%     
    for k=1:N
        Cost_function=Cost_function+0.5*(x{k})'*Q*x{k}+0.5*(u{k})'*R*u{k};
        x{k+1}=RHO*x{k}+Nu*u{k}+Ne*e;
        y{k+1}=Cy*x{k}+Du*u{k};
%         constrain=[constrain,x{k+1}==kest.A*x{k}+kest.B*u{k}+Ne*((((Cy*RHO*Cy')/(Cy*Cy'))^(k-1))*e+v_cumu)];
%         constrain=[constrain,x{N+1}==0];
%         constrain=[constrain,x{N+1}<=0.0008];
%         constrain=[constrain,x{N+1}>=-0.0008];
%         constrain=[constrain,x{k+1}==RHO*x{k}+Nu*u{k}+Ne*e];
%         constrain=[constrain,y{k+1}==Cy*x{k}+Du*u{k}];
%         v_cumu=(Cy*RHO*Cy')/(Cy*Cy')*v_cumu+v(:,k);
    end
    optimize(constrain,Cost_function,sdpsettings('solver','quadprog'));
    if j==1
        u_prev=0;
    else
        u_prev=u_in;
    end
    
    u_in=double(u{1});
    y_hat=Cy*double(x{2})+Du*u_in;
    input=[input,u_in];

    [X,Y,Z]=real_system(X,u_in,RHO,Nu,Nd,Cz,Dd,Cy,Du,ground(j),v2);

    x0=X;
    
    % Computational time
    diff_time=cputime-t;
    com_Time=[com_Time,diff_time];
    
    % Compute the error
    e=Y-Cy*x0-Du*u_in;
    
    % Lists for plotting
    energy=[energy,-1/2*L*u_in'*(u_in-u_prev)/dt];
    error=[error,e];
    state(1:(2*cas),j+1)=X;
    accel(1:(cas),j+1)=Y;
end

%% Non linear part 
input=[];
init=["e",0,0,0,0]
for k=1:it
    k
    if k>=2
        u_in=U(1);
        [X,Y,Z]=NL_system(init,X,u_in,RHO,RHO2,Nu,Nu2,Nd,Nd2,Dd,Cy,Cz,Du,ground(k),v2,Ne1,Ne2,0,P12,P1,1);
        init=[ed,x_out_of_deform,x_prev,yy,x_in]
        %  Compute the error
        e=Y-Cy*x0-Du*u_in;
        state(1:(2*cas),k+1)=X;
        input=[input;u_in];
        % Lists for plotting
%         energy=[energy,-1/2*L*u_in'*(u_in-u_prev)/dt];
        error=[error,e];
    end
    
    U=fmincon(@(u)cost(init,X,u,RHO,RHO2,Nu,Nu2,Nd,Nd2,Dd,Cy,Cz,Du,0,0,Ne1,Ne2,e,P12,P1,Q,R,N),zeros(1,N),[],[],[],[],[],[])%,@(u)constraints(init,X1,u,RHO,RHO2,Nu,Nu2,Nd,Nd2,Dd,Cy,Cz,Du,0,v,Ne1,Ne2,e,P12,P1,N));
end

%%
v=1/769*randn(cas,it+1);
v_cumu=v(:,1);

accel=[0];
e=0;
Y=0;
X1=zeros(2*cas,1);
X=X1;
init=["e",0,0,0,0]
for j=1:it
%     if j==1
%         Vars=whos;
%         PersistentVars=Vars([Vars.global]);
%         PersistentVarNames={PersistentVars.name};
%         clear(PersistentVarNames{:});
%     end
%     j;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEST OF THE KALMAN FILTER

    X=RHO*X+Ne*(e)
    Y=Cy*X
    v_cumu=(Cy*RHO*Cy')/(Cy*Cy')*v_cumu+v(j+1);
    state2(1:(2*cas),j+1)=X; %To plot the Falman filter result 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SIMULATIONS OF SYSTEMS WITHOUT CONTROLLER

%NOISE
%     v2=10/769*randn(cas,1);

    % SYSTEM WITHOUT PLASTICITY
%     [X_without,Y_without,Z_without]=real_system(X1,0,RHO,Nu,Nd,Cz,Dd,Cy,Du,ground(j),v2);

    % SYSTEM WITH PLASTICITY
%     [X_without,Y_without,Z_without]=NL_system(init,X1,0,RHO,RHO2,Nu,Nu2,Nd,Nd2,Dd,Cy,Cz,Du,ground(j),v2,Ne1,Ne2,0,P12,P1,j);
%     init=[ed,x_out_of_deform,x_prev,yy,x_in]  

%     X1=X_without;

%     accel=[accel,Y_without];

%     e=(Y_without-Y); 
% error between the kalman filter and the real system

% list for plotting the system without controller
%     normal_state(1:(2*cas),j+1)=X_without;
%     accel_without(1:(cas),j+1)=Y_without;
%     
% %     e=Y_without-(Cy*X+v(j+1))
   
end

%% PLOT PART

t=(0:it)*dt;
figure
% subplot(4,1,1)
plot(state(1,:),'Color','k')
hold on
plot(normal_state(1,:),'Color','r')
title("Displacement of the first floor with non linearity")%with and without controller")

% IN CASE OF A SECOND FLOOR
% subplot(4,1,2)
% plot(t,state(2,:),'Color','k')
% hold on
% plot(t,normal_state(2,:),'Color','r')
% title("Celerity of the first floor with and without controller")

%INPUT
% subplot(4,1,3)
% plot(input(1,:))
% title("Input (displacement of the mass of the tendon system for the first floor)")
% subplot(5,1,4)
% plot(input(2,:))
% title("Input (displacement of the mass of the tendon system for the second floor)")
% subplot(4,1,4)
% plot(t(2:size(t,2)),energy(1,:))
% title("Power necessary to apply the input")

%INPUT IN BIGGER WINDOW
figure 
plot(input)
% 
%ERROR
% % figure
% % plot(t,error)
% 
%  figure
%  plot(t,accel,'Color','k')
% hold on
% plot(t,accel_without,'Color','r')
% 
% figure 
% plot(t(2:it+1),com_Time)
% 
%figure
% plot(state2(1,:))


