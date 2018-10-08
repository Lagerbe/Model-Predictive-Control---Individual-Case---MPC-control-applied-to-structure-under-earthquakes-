function [k,y,y2]=hysteresis(init,x,ke,kd,Re,K1)
    global ed x_out_of_deform x_prev yy x_in
%% Initialization
    if isempty(x_out_of_deform)
        ed="e";
        x_prev=0;
        x_out_of_deform=0;
        yy=0;
        x_in=0;
    end
    if ~isempty(init)
        ed=init(1);
        x_prev=str2double(init(3));
        x_out_of_deform=str2double(init(2));
        yy=str2double(init(4));
        x_in=str2double(init(5));
    end
    %% Elastic deformation
    if ed=="e" & x_out_of_deform==0 & abs(ke*(x-x_out_of_deform))<=Re 
        "cas1 elastic 1"
        k=K1;
        x_prev=x;
        y=yy;
    elseif ed=="e" & x_out_of_deform~=0 & abs(ke*(x-x_out_of_deform))<=2*Re 
        "cas2 elastic"
        k=K1;
        x_prev=x;
        y=yy;
 %% Elastic to plastic deformation
    elseif ed=="e" & x_out_of_deform==0 & abs(ke*(x-x_out_of_deform))>Re
        "cas3 plastic 1"
        k=kd;
        ed="d";
        x_prev=x;
        x_in=x_prev;
        yy=yy+(K1-kd)*x_in;%K1*(x-x_out_of_deform)
        y=yy;
    elseif ed=="e" & x_out_of_deform~=0 & abs(ke*(x-x_out_of_deform))>2*Re
        "cas4 plastic"
        k=kd;
        ed="d";
        x_prev=x;
        x_in=x_prev;
        yy=yy+(K1-kd)*x_in;
        y=yy;
 %% Evolution in the deformation
 % Extention
    elseif ed=="d" & x_out_of_deform==0 & ke*(x_prev-x_out_of_deform)>Re & x-x_prev>=0
        "cas5 elastic to plastic 1"
        k=kd;
        ed="d";
        x_prev=x;
        y=yy;
 %Compression
    elseif ed=="d" & x_out_of_deform==0 & ke*(x_prev-x_out_of_deform)<-Re & x-x_prev<=0
        "cas6 elastic to plastic"
        k=kd;
        ed="d";
        x_prev=x;
        y=yy;
 % Extention
    elseif ed=="d" & x_out_of_deform~=0 & ke*(x_prev-x_out_of_deform)>2*Re & x-x_prev>=0
        "cas7 elastic to plastic"
        k=kd;
        ed="d";
        x_prev=x;
        y=yy;
 %Compression
    elseif ed=="d" & x_out_of_deform~=0 & ke*(x_prev-x_out_of_deform)<-2*Re & x-x_prev<=0
        "cas8 elastic to plastic"
        k=kd;
        ed="d";
        x_prev=x;
        y=yy;
   %% Plastic to elastic deformation
    %Extention to elasticity
    elseif ed=="d" & x_out_of_deform==0 & ke*(x_prev-x_out_of_deform)>Re & x-x_prev<0 & abs(ke*(x-x_prev))<2*Re
        "cas9 plastic extention to elastic 1"
        k=K1;
        ed="e";
        x_out_of_deform=x_prev;
        x_prev=x;
        yy=yy+(kd-K1)*x_out_of_deform;
        y=yy;
    %Compression to elasticity
   elseif ed=="d" & x_out_of_deform==0 & ke*(x_prev-x_out_of_deform)<-Re & x-x_prev>0 & abs(ke*(x-x_prev))<2*Re
       "cas10 plastic compression to elastic 1"
        k=K1;
        ed="e";
        x_out_of_deform=x_prev;
        x_prev=x;
        %kd*(x-x_in)+K1*x_out_of_deform
        yy=yy+kd*(x-x_in)+K1*x_out_of_deform
        %(kd-K1)*x_out_of_deform;
        y=yy;
   %Extention to elasticity
   elseif ed=="d" & x_out_of_deform~=0 & ke*(x_prev-x_out_of_deform)>2*Re & x-x_prev<0 & abs(ke*(x-x_prev))<2*Re
       "cas11 plastic extention to elastic"
        k=K1;
        ed="e";
        x_out_of_deform=x_prev;
        x_prev=x;
        yy=yy+(kd-K1)*x_out_of_deform;
        y=yy;
   %Compression to elasticity
   elseif ed=="d" & x_out_of_deform~=0 & ke*(x_prev-x_out_of_deform)<-2*Re & x-x_prev>0 & abs(ke*(x-x_prev))<2*Re
       "cas12 plastic compression to elastic"
        k=K1;
        ed="e";
        x_out_of_deform=x_prev;
        x_prev=x;
        yy=yy+(kd-K1)*x_out_of_deform;
        y=yy;
        %% IN CASE OF A SMALL SAMPLING TIME
   elseif ed=="d" & x_out_of_deform==0 & ke*(x_prev-x_out_of_deform)>Re & x-x_prev<0 & abs(ke*(x-x_prev))>2*Re
        k=K1;
        ed="e";
        x_out_of_deform=x_prev;
        x_prev=x;
        yy=yy+(kd-K1)*x_out_of_deform-(kd-K1)*x;
        y=yy;
    %Compression to elasticity
   elseif ed=="d" & x_out_of_deform==0 & ke*(x_prev-x_out_of_deform)<-Re & x-x_prev>0 & abs(ke*(x-x_prev))>2*Re
        k=K1;
        ed="e";
        x_out_of_deform=x_prev;
        x_prev=x;
        yy=yy+(kd-K1)*x_out_of_deform-(kd-K1)*x;
        y=yy;
   %Extention to elasticity
   elseif ed=="d" & x_out_of_deform~=0 & ke*(x_prev-x_out_of_deform)>2*Re & x-x_prev<0 & abs(ke*(x-x_prev))>2*Re
        k=K1;
        ed="e";
        x_out_of_deform=x_prev;
        x_prev=x;
        yy=yy+(kd-K1)*x_out_of_deform-(kd-K1)*x;
        y=yy;
   %Compression to elasticity
   elseif ed=="d" & x_out_of_deform~=0 & ke*(x_prev-x_out_of_deform)<-2*Re & x-x_prev>0 & abs(ke*(x-x_prev))>2*Re
        k=K1;
        ed="e";
        x_out_of_deform=x_prev;
        x_prev=x;
        yy=yy+(kd-K1)*x_out_of_deform-(kd-K1)*x;
        y=yy;
%     else
%         k=K1;
%         y=yy;
%         x_prev=x;
    end
    y2=y+k*x;
%     evalin('base','yy=yy;')