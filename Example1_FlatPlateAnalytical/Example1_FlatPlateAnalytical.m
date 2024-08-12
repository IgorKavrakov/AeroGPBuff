function [FP] = Example1_FlatPlateAnalytical(FP,Lag)
% This function gives the training and prediction data based on analytical flat plate solution

% By Igor Kavrakov

%%%%%%%%% COPYRIGHT NOTICE %%%%%%%%% 
%  This file is part of AeroGPBuff.
%  AeroGPBuff  is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation, either version 3 of the License, or
%  (at your option) any later version.
% 
%  AeroGPBuff  is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
% 
%  You should have received a copy of the GNU General Public License
%  along with AeroGPBuff .  If not, see <https://www.gnu.org/licenses/>.

% Copyright (c) Igor Kavrakov, Guido Morgenthal, Allan McRobie 2024

%% Training
ds=FP.Par.ds;
if ~isfield(FP,'Pred')||strcmp(FP.Pred.Excitation,'Train') %Training and prediction set are same
VrR=FP.Train.Vr;
Alpha_AmpR=FP.Train.Alpha_AmpRel;

Cycl=FP.Train.Cycl;
TPoint=Cycl*VrR(end)./ds;
tau=ds*(0:1:TPoint-1).'; %Reduced time

Vr=ds*TPoint./(1:TPoint)';
Vr(TPoint/2:end)=-inf;
Vals= Vr>=VrR(1) & Vr<=VrR(2);
    alpha_h=[];
    alpha_h_d=[];
    alpha_a=[];
    alpha_a_d=[];
    CL=[];
    CM=[];

for i=1:length(FP.Train.Alpha_STD)
    Alpha_STD=FP.Train.Alpha_STD(i)*pi/180;

    X_FFT=(1-Alpha_AmpR).*rand(TPoint,1) + Alpha_AmpR;
    Phi_rand = exp(1j.*2*pi*rand(TPoint,1));
    X_FFT(~Vals)=0;

    alpha_h_l=real(ifft(X_FFT.*Phi_rand));
    alpha_h_d_l=real(ifft(1j*2*pi./Vr.*X_FFT.*Phi_rand));

    X_FFT=(1-Alpha_AmpR).*rand(TPoint,1) + Alpha_AmpR;
    Phi_rand = exp(1j.*2*pi*rand(TPoint,1));
    X_FFT(~Vals)=0;

    alpha_a_l=real(ifft(X_FFT.*Phi_rand));
    alpha_a_d_l=real(ifft(1j*2*pi./Vr.*X_FFT.*Phi_rand));
    alpha_a_2d=real(ifft((1j*2*pi./Vr).^2.*X_FFT.*Phi_rand));


    FactAmp_h=Alpha_STD./std(alpha_h_l); %This factor is close -> but not exactly the same as the input, but it does not matter!
    FactAmp_a=Alpha_STD./std(alpha_a_l);

    alpha_h_l=(alpha_h_l-mean(alpha_h_l)).*FactAmp_h; 
    alpha_h_d_l=(alpha_h_d_l-mean(alpha_h_d_l)).*FactAmp_h; 

    alpha_a_l=(alpha_a_l-mean(alpha_a_l)).*FactAmp_a; 
    alpha_a_d_l=(alpha_a_d_l-mean(alpha_a_d_l)).*FactAmp_a; 
    alpha_a_2d=(alpha_a_2d-mean(alpha_a_2d)).*FactAmp_a; 

    Phi_se=1-0.165.*exp(-0.089.*tau)-0.335.*exp(-0.6.*tau);% Wagner function        
    h_primeprimeB=sec(alpha_h_l).^2.*alpha_h_d_l;     %H_dd/B. Relationship 1+tan^2(alpha)=sec^2(alpha) is used.
    h_conv=ConvFFT(Phi_se,h_primeprimeB,ds);      %Rise due to vert vel
    a_conv=ConvFFT(Phi_se,alpha_a_d_l,ds);      %Rise due to rotation 
    a_d_conv=ConvFFT(Phi_se,alpha_a_2d,ds/4); %Rise due to angular velocity

    CLh=2*pi*(-h_conv-1/4*h_primeprimeB);                                      %Lift due to h
    CMh=pi/2*( h_conv);                                                       %Moment due to h
    CLa=2*pi*(-a_conv-a_d_conv-1/4*alpha_a_d_l);                   %Lift due to a
    CMa=pi/2*( a_conv+a_d_conv-1/4*alpha_a_d_l-alpha_a_2d/32);  %Moment due to a

    alpha_h_l(:,2:Lag)=0; 
    alpha_h_l(:,2:Lag)=0; 

    for ll=1:Lag %Regressors
        alpha_h_l(1+ll:end,1+ll)=alpha_h_l(1:TPoint-ll,1);   
        alpha_a_l(1+ll:end,1+ll)=alpha_a_l(1:TPoint-ll,1); 
    end
    
    alpha_h   = [alpha_h;atan(alpha_h_l)];
    alpha_h_d = [alpha_h_d;alpha_h_d_l];      % Concattanate
    alpha_a   = [alpha_a;alpha_a_l];
    alpha_a_d = [alpha_a_d;alpha_a_d_l];      % Concattanate
 
    CL= [CL;CLh-mean(CLh)+CLa-mean(CLa)];
    CM= [CM;CMh-mean(CMh)+CMa-mean(CMa)];

    
end
    TPoint= TPoint*length(FP.Train.Alpha_STD);
    tau=ds*(0:1:TPoint-1).'; %Reduced time
    
    % Buffeting training
    if isfield(FP.Train,'Buf')&&FP.Train.Buf==1

        r    = FP.Par.r;
        Iw   = FP.Par.Iw;
        Tw   = FP.Par.Tw;    
        Cyclb= FP.Train.Cycl_b; 
        Vr_b = FP.Train.Vr_b;
        tau_max = Vr_b * Cyclb;

        if ~isfield(FP.Train,'Coupled')||FP.Train.Coupled~=1
        tau_max = max(tau_max,r/0.01);
        else
        TPointS=Cycl*VrR(end)./ds;
        tau_max = TPointS*ds;
        end

        tau_b = 0:ds:tau_max-ds;
        TPointB = length(tau_b);    
        dk= 1/tau_max;
        kw=(1:1:TPointB/2).*dk;  
        tau_b=ds*(0:1:TPointB-1).'; %Reduced time

        alpha_w=[];
        alpha_w_d=[];
        CLw=[];
        CMw=[];
        Sw=[];
        for i=1:length(Iw)
            [~,Sw1,w,w_conv] = FP_Analy_Buf (Iw(i),Tw,kw,r,TPointB);
            Sw = [Sw Sw1'];     
            alpha_w_l = atan(w');
            alpha_w_l(:,2:Lag)=0; 
            for ll=1:Lag %Regressors
                alpha_w_l(1+ll:end,1+ll)=alpha_w_l(1:TPointB-ll,1);   
            end

            alpha_w_d_l=(alpha_w_l(2:end,1)-alpha_w_l(1:end-1,1))./ds;
            alpha_w_d_l(end+1)=0;
            alpha_w_d = [alpha_w_d;alpha_w_d_l];      % Concattanate
            alpha_w   = [alpha_w;alpha_w_l];
            CLw= [CLw;-2*pi*(w_conv-mean(w_conv))];
            CMw= [CMw;pi/2*(w_conv-mean(w_conv))];

        end
            TPointB = TPointB*length(Iw);
            tau_b=ds*(0:1:TPointB-1).'; %Reduced time
    else
        alpha_w = [];
        alpha_w_d = [];
        CLw = [];
        CMw = [];
        TPointB = [];
        tau_b = [];
        Sw = [];
        kw = [];

    end
%Training
FP.Train.alpha_h=alpha_h;
FP.Train.alpha_h_d=alpha_h_d;
FP.Train.alpha_a=alpha_a;
FP.Train.alpha_a_d=alpha_a_d;
FP.Train.CL=CL;
FP.Train.CM=CM;
FP.Train.Samp=TPoint;
FP.Train.tau=tau;

FP.Train.alpha_w=alpha_w;
FP.Train.alpha_w_d=alpha_w_d;
FP.Train.CLw=CLw;
FP.Train.CMw=CMw;
FP.Train.Samp_b=TPointB;
FP.Train.tau_b=tau_b;
FP.Train.Sw = Sw;
FP.Train.k  = kw;    
    
end

%% Prediction 
    if ~isfield(FP,'Pred')||~isfield(FP.Pred,'Excitation')||strcmp(FP.Pred.Excitation,'Train') %Training and prediction set are same
    alpha_h_pr=alpha_h;
    alpha_h_d_pr=alpha_h_d;
    alpha_a_pr=alpha_a;
    alpha_a_d_pr=alpha_a_d;
    CL_pr=CL;
    CM_pr=CM;
    PPoint=TPoint;
    
    alpha_w_pr = alpha_w;
    alpha_w_d_pr = alpha_w_d;
    CLw_pr = CLw;
    CMw_pr = CMw;
    PPointB = TPointB;
    tau_b_pr = tau_b; 
      
    elseif strcmp(FP.Pred.Excitation,'Rand')%Random excitation
        load('FP_RandInp.mat','a_pr','a_d_pr','a_2d_pr','h_pr','h_d_pr','ds'); %Loads the input for the prediction ,alpha_pr, alpha_d_pr, and alpha_2d_pr, ds(Reduced time step). These are used for prediction    
        PPoint=length(a_pr); 
        tau_predict=(0:ds:ds*(PPoint-1))'; %Time & Reduced time
        Phi_se=1-0.165.*exp(-0.089.*tau_predict)-0.335.*exp(-0.6.*tau_predict);% Wagner function
        h_conv=ConvFFT(Phi_se,h_d_pr,ds);      %Rise due to vert vel
        a_conv=ConvFFT(Phi_se,a_d_pr,ds);  %Rise lift due to rotation 
        a_d_conv=ConvFFT(Phi_se,a_2d_pr,ds/4); %Rise lift due to angular velocity

        alpha_h_pr=zeros(PPoint,1+Lag); 
        alpha_h_d_pr=1./(1+h_pr.^2).*h_d_pr;
        alpha_a_pr=zeros(PPoint,1+Lag); 
        alpha_a_d_pr=a_d_pr;
        
        CMh_pr=pi/2*( h_conv);                                                     
        CMa_pr=pi/2*( a_conv+a_d_conv-1/4*a_d_pr-a_2d_pr/32);
              
        CL_pr=zeros(PPoint,1);
        CM_pr=CMh_pr+CMa_pr;

        for ll=0:Lag% Add the lag, leaving out 0 initial conditions for the regressors
            alpha_h_pr(1+ll:end,1+ll)  =atan(h_pr(1:end-ll,1));
            alpha_a_pr(1+ll:end,1+ll)  =a_pr(1:end-ll,1);              
        end
        
    elseif strcmp(FP.Pred.Excitation,'SinW')
        Vr=FP.Pred.VrW;
        Alpha_Amp=FP.Pred.Alpha_AmpW*pi/180;
        Cycl=FP.Pred.CyclW;  
        
        PPointB=sum(Cycl.*Vr./ds+1)*length(Alpha_Amp);
        
        alpha_w_pr=zeros(PPointB,1+Lag); 
        alpha_w_d_pr=zeros(PPointB,1);

        CLw_pr=zeros(PPointB,1); 
        CMw_pr=zeros(PPointB,1); 
        
        Per=0;Num=0;Inc=sum((Cycl.*Vr./ds+1)); %Initiate
        for i=1:length(Alpha_Amp)
            for j=1:length(Vr)
              
              K=2*pi/Vr(j); k=K/2; %Reduced frequency & Halfreduced frequency
              tau_b_pr=(0:ds:(Vr(j)*Cycl(j))).'; %Total reduced time for number of Cycl
              alpha=Alpha_Amp(i)*exp(1j.*K.*tau_b_pr-1j*pi/2); %Sinusoidal alpha          

              NElemVR=length(alpha);

              Sears=0.5*0.13^2./(0.13^2+k.^2)+0.5./(1+k.^2)-1j.*k.*(0.5*0.13./(0.13^2+k.^2)+0.5./(1+k.^2));                    
              CLw_pr(1+Num+Per:Per+Num+NElemVR)=-2*pi*real(Sears*alpha);       %Lift coefficient due to vert gust
              CMw_pr(1+Num+Per:Per+Num+NElemVR)=pi/2*real(Sears*alpha);       %Moment coefficient due to vert gust

              alpha=real(alpha);
              NElemVR=length(alpha);

              alpha_w_d_pr(1+Num+Per:Per+Num+NElemVR-1,1) = (alpha(2:end) - alpha(1:end-1))/ds;
                  for ll=0:Lag %Regressors
                      alpha_w_pr(1+ll+Num+Per:Per+Num+NElemVR,1+ll)   =atan(alpha(1:NElemVR-ll,1));
                      alpha_w_pr(1+Num+Per:1+ll+Num+Per,1+ll)         =atan(alpha(NElemVR-ll:end,1));
                  end
              Num=Num+NElemVR; %Reset               
            end 
            Per=Per+Inc;Num=0; %Reset
        end
        alpha_h_pr=[];  
        alpha_h_d_pr=[];
        alpha_a_pr=[]; 
        alpha_a_d_pr=[]; 
        CL_pr=[];
        CM_pr=[];
        PPoint=[];        

        
    elseif strcmp(FP.Pred.Excitation,'SinH')|| strcmp(FP.Pred.Excitation,'SinA')  %Sinusoidal excitation
        Vr=FP.Pred.Vr;
        Alpha_Amp=FP.Pred.Alpha_Amp*pi/180;
        Cycl=FP.Pred.Cycl;

        PPoint=sum(Cycl.*Vr./ds+1)*length(Alpha_Amp);
        
        alpha_h_pr=zeros(PPoint,1+Lag); 
        alpha_h_d_pr=zeros(PPoint,1);

        alpha_a_pr=zeros(PPoint,1+Lag); 
        alpha_a_d_pr=zeros(PPoint,1);
        
        CLh_pr=zeros(PPoint,1); 
        CMh_pr=zeros(PPoint,1); 
        CLa_pr=zeros(PPoint,1); 
        CMa_pr=zeros(PPoint,1); 

         Per=0;Num=0;Inc=sum((Cycl.*Vr./ds+1)); %Initiate
        for i=1:length(Alpha_Amp)
            for j=1:length(Vr)

            [alphaVr,alphaVr_d,CLhVr,CMhVr,CLaVr,CMaVr] = FP_Analy_SE (Vr(j),ds,Cycl(j),Alpha_Amp(i));
            NElemVR=length(alphaVr);
            
            if strcmp(FP.Pred.Excitation,'SinH')
            alpha_h_d_pr(1+Num+Per:Per+Num+NElemVR,1) =1./(1+alphaVr.^2).*alphaVr_d;
            CLh_pr(1+Num+Per:Per+Num+NElemVR)=CLhVr;%Lift coefficient due to vert disp
            CMh_pr(1+Num+Per:Per+Num+NElemVR)=CMhVr;%Moment coefficient due to rot 
            elseif strcmp(FP.Pred.Excitation,'SinA')     
            alpha_a_d_pr(1+Num+Per:Per+Num+NElemVR,1) =alphaVr_d;               
            CLa_pr(1+Num+Per:Per+Num+NElemVR)=CLaVr;%Lift coefficient due to vert disp
            CMa_pr(1+Num+Per:Per+Num+NElemVR)=CMaVr;%Moment coefficient due to rot
            end
              for ll=0:Lag %Regressors
                  if strcmp(FP.Pred.Excitation,'SinH')
                  alpha_h_pr(1+ll+Num+Per:Per+Num+NElemVR,1+ll)   =atan(alphaVr(1:NElemVR-ll,1));
                  alpha_h_pr(1+Num+Per:1+ll+Num+Per,1+ll)         =atan(alphaVr(NElemVR-ll:end,1));
                  elseif strcmp(FP.Pred.Excitation,'SinA')
                  alpha_a_pr(1+ll+Num+Per:Per+Num+NElemVR,1+ll)   =alphaVr(1:NElemVR-ll,1);
                  alpha_a_pr(1+Num+Per:1+ll+Num+Per,1+ll)         =alphaVr(NElemVR-ll:end,1);
                  end
              end
              Num=Num+NElemVR; %Reset
            end 
            Per=Per+Inc;Num=0; %Reset
        end
      
        CL_pr=CLh_pr+CLa_pr;
        CM_pr=CMh_pr+CMa_pr;
        
        alpha_w_pr=alpha_h_pr*0;  
        alpha_w_d_pr=alpha_h_d_pr*0; 
        CLw_pr=CL_pr*0;
        CMw_pr=CM_pr*0;
        PPointB=PPoint;
        tau_b_pr=[];
        
    elseif strcmp(FP.Pred.Excitation,'Stat') %Static excitation
        Alpha_Amp=FP.Pred.StatAlpha*pi/180;
        PPoint=(FP.Pred.StatTau./ds+1)*length(Alpha_Amp);
        
        alpha_h_pr=zeros(PPoint,1+Lag); 
        alpha_h_d_pr=zeros(PPoint,1);

        alpha_a_pr=zeros(PPoint,1+Lag); 
        alpha_a_d_pr=zeros(PPoint,1);
        
        CL_pr=zeros(PPoint,1);
        CM_pr=zeros(PPoint,1);    
        
        Inc=(FP.Pred.StatTau./ds+1); %Initiate
        for i=1:length(Alpha_Amp)
           alpha_a_pr(1+(i-1)*Inc:i*Inc,:)=Alpha_Amp(i);
           CL_pr(1+(i-1)*Inc:i*Inc)=Alpha_Amp(i)*2*pi; %Slope 2*pi
           CM_pr(1+(i-1)*Inc:i*Inc)=Alpha_Amp(i)*pi/2; %Slope pi/2            
        end
    
        alpha_w_pr=alpha_h_pr*0;  
        alpha_w_d_pr=alpha_h_d_pr*0; 
        CLw_pr=CL_pr*0;
        CMw_pr=CM_pr*0;
        PPointB=PPoint;
        tau_b_pr=FP.Pred.StatTau;
        
    elseif strcmp(FP.Pred.Excitation,'Aeroelastic')  %Flutter       
        M=[FP.Par.m FP.Par.I];
        C=4*pi*[FP.Par.fh FP.Par.fa].*FP.Par.psi.*M;
        K=4*pi^2*[FP.Par.fh FP.Par.fa].^2.*M;

        M=M.*FP.Par.U^2/FP.Par.B^2; % Analysis w.r.t. nondimensional time
        C=C.*FP.Par.U/FP.Par.B;     % Analysis w.r.t. nondlimensional time
        
        PPoint=floor(FP.Par.Redtime./ds);
        tau_predict=(0:ds:ds*(PPoint-1))'; %Time & Reduced time
        Phi_se=1-0.165.*exp(-0.089.*tau_predict)-0.335.*exp(-0.6.*tau_predict);% Wagner function
        
        u=zeros(PPoint,2); %Dispalcements
        u_d=zeros(PPoint,2); %Velocity
        u_2d=zeros(PPoint,2); %Acceleration
        
        u(1,1)=0.5; %Initial conditions
        
        CL_pr=zeros(PPoint,1);
        CM_pr=zeros(PPoint,1);
        Pres=1/2*FP.Par.rho*FP.Par.U^2*FP.Par.B.*[1 FP.Par.B];
        B=FP.Par.B;
        
       for i=1:PPoint %Time integration loop
        [~,h_conv]=ConvFFT(Phi_se(1:i),u_2d(1:i,1)./B,ds);
        [~,a_conv]=ConvFFT(Phi_se(1:i),u_d(1:i,2),ds);
        [~,a_d_conv]=ConvFFT(Phi_se(1:i),u_2d(1:i,2).*1/4,ds);
        Conv=h_conv+a_conv+a_d_conv;
                
        CL_pr(i)=2*pi*(-1/4*u_2d(i,1)/B-1/4*u_d(i,2)-Conv);
        CM_pr(i)=pi/2*(-1/4*u_d(i,2)-u_2d(i,2)/32   +Conv);

        p=Pres.*[CL_pr(i) CM_pr(i)];
        [u(i+1,:),u_d(i+1,:),u_2d(i+1,:)] = NewmarkSDOF(p,u(i,:),u_d(i,:),u_2d(i,:),1/4,1/2,ds,M,C,K);    
       end
       %Store extra
        FP.Pred.u=u; FP.Pred.u_d=u_d; FP.Pred.u_2d=u_2d;
        FP.Par.K=K; FP.Par.C=C;       FP.Par.M=M;
        FP.Par.Pres=Pres;
        
        alpha_h_pr=[];  
        alpha_h_d_pr=[];
        alpha_a_pr=[]; 
        alpha_a_d_pr=[];
        
        alpha_w_pr=[];
        alpha_w_d_pr=[];
        CLw_pr=[];
        CMw_pr=[];
        PPointB=[];
        tau_b_pr=[];
        
    elseif strcmp(FP.Pred.Excitation,'Buffeting') %Buffeting
        M=[FP.Par.m FP.Par.I];
        C=4*pi*[FP.Par.fh FP.Par.fa].*FP.Par.psi.*M;
        K=4*pi^2*[FP.Par.fh FP.Par.fa].^2.*M;

        M=M.*FP.Par.U^2/FP.Par.B^2; % Analysis w.r.t. nondimensional time
        C=C.*FP.Par.U/FP.Par.B;     % Analysis w.r.t. nondlimensional time

        %Parameters & time
        PPoint=floor(FP.Par.Redtime./ds);
        tau_predict=(0:ds:ds*(PPoint-1))'; %Time & Reduced time
        PPointB=PPoint;
        tau_b_pr=tau_predict;
    
        % Wagner function - self  excited
        Phi_se=1-0.165.*exp(-0.089.*tau_predict)-0.335.*exp(-0.6.*tau_predict);% Wagner function
       
        %Buffeting
        r    = FP.Par.r;
        Iw   = FP.Par.Iw;
        Tw   = FP.Par.Tw;
        PPoint=floor(FP.Par.Redtime./ds);
        dk= 1/FP.Par.Redtime;
        kw=(1:1:PPoint/2).*dk;
        [~,~,w,w_conv] = FP_Analy_Buf(Iw,Tw,kw,r,PPoint); % w is w/U
       
        %Buffeting LU 
        u=zeros(PPoint,2); %Dispalcements
        u_d=zeros(PPoint,2); %Velocity
        u_2d=zeros(PPoint,2); %Acceleration
        
        CL_pr=zeros(PPoint,1);
        CM_pr=zeros(PPoint,1);
        Pres=1/2*FP.Par.rho*FP.Par.U^2*FP.Par.B.*[1 FP.Par.B];
        B=FP.Par.B;
        
        alpha_w_pr = atan(w');
        alpha_w_pr(:,2:Lag+1)=0; 
        for ll=1:Lag %Regressors
            alpha_w_pr(1+ll:end,1+ll)=alpha_w_pr(1:PPointB-ll,1);   
        end                      
        alpha_w_d_pr=(alpha_w_pr(2:end,1)-alpha_w_pr(1:end-1,1))./ds;
        alpha_w_d_pr(end+1)=0;        
        
        CLw_pr=-2*pi*w_conv;
        CMw_pr=pi/2*w_conv;

       for i=1:PPoint %Time integration loop
        [~,h_conv]=ConvFFT(Phi_se(1:i),u_2d(1:i,1)./B,ds);
        [~,a_conv]=ConvFFT(Phi_se(1:i),u_d(1:i,2),ds);
        [~,a_d_conv]=ConvFFT(Phi_se(1:i),u_2d(1:i,2).*1/4,ds);
        Conv=h_conv+a_conv+a_d_conv;
                
        CL_pr(i)=2*pi*(-1/4*u_2d(i,1)/B-1/4*u_d(i,2)-Conv);
        CM_pr(i)=pi/2*(-1/4*u_d(i,2)-u_2d(i,2)/32   +Conv);

        p=Pres.*[CL_pr(i) CM_pr(i)] + Pres.*[CLw_pr(i) CMw_pr(i)];
        [u(i+1,:),u_d(i+1,:),u_2d(i+1,:)] = NewmarkSDOF(p,u(i,:),u_d(i,:),u_2d(i,:),1/4,1/2,ds,M,C,K);    
       end
       
       %Store extra
        FP.Pred.u=u;       FP.Pred.u_d=u_d;       FP.Pred.u_2d=u_2d;
        FP.Par.K=K; FP.Par.C=C;       FP.Par.M=M;
        FP.Par.Pres=Pres;

        alpha_h_pr=[];  
        alpha_h_d_pr=[];
        alpha_a_pr=[]; 
        alpha_a_d_pr=[];

        
    end

%% Store
%Prediction
FP.Pred.alpha_h=alpha_h_pr;
FP.Pred.alpha_h_d=alpha_h_d_pr;
FP.Pred.alpha_a=alpha_a_pr;
FP.Pred.alpha_a_d=alpha_a_d_pr;
FP.Pred.CL=CL_pr;
FP.Pred.CM=CM_pr;
FP.Pred.Samp=PPoint;


FP.Pred.alpha_w=alpha_w_pr;
FP.Pred.alpha_w_d=alpha_w_d_pr;
FP.Pred.CLw=CLw_pr;
FP.Pred.CMw=CMw_pr;
FP.Pred.Samp_b=PPointB;
FP.Pred.tau_b=tau_b_pr;

    
end

%% Self-excited forces due to sinusoidal effective angle
function [alpha,alpha_d,CLh,CMh,CLa,CMa] = FP_Analy_SE (Vr,ds,Cycl,Amp)
          K=2*pi/Vr; k=K/2; %Reduced frequency & Halfreduced frequency
          tau=(0:ds:(Vr*Cycl)).'; %Total reduced time for number of Cycl
          alpha=Amp*exp(1j.*K.*tau-1j*pi/2); %Sinusoidal alpha          
          alpha_d=(1j.*K.*Amp*exp(1j.*K.*tau-1j*pi/2));     %Angular velocity 
          alpha_2d=((1j.*K).^2.*Amp*exp(1j.*K.*tau-1j*pi/2));%Angular acceleration 
          
          Theo=conj(besselh(1,k))./(conj(besselh(1,k))+ 1j.*conj(besselh(0,k)));%Theodorsen function
          
          CLh=2*pi*real(-Theo*alpha-1/4*alpha_d);                     %Lift coefficient due to vert disp
          CMh=pi/2*real( Theo*alpha);                                      %Moment coefficient due to rot  
          CLa=2*pi*real(-Theo*(alpha+alpha_d*1/4)-1/4*alpha_d);            %Lift coefficient due to vert disp
          CMa=pi/2*real( Theo*(alpha+alpha_d*1/4)-1/4*alpha_d-alpha_2d/32);%Moment coefficient due to rot         
         
          alpha=real(alpha);
          alpha_d=real(alpha_d);
end

%% Fast convolution
function [Conv,Cent] = ConvFFT(M1,M2,dt)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
M1 = [M1;zeros(length(M1),1)];
M2 = [M2;zeros(length(M2),1)];
  
Conv = ifft(fft(M1,[],1).*fft(M2,[],1),[],1);
Conv = Conv(1:length(M1)/2,:)*dt;
Cent=Conv(end,:);
end

%% Buffeting forces - von Karman
function [Sww,Sw,w,w_conv,Sears] = FP_Analy_Buf (Iw,Tw,kw,r,TPointB)
        Sww=4.*Iw.^2 .* Tw .* (1+755.2.*(kw.*r).^2)./(1+283.2.*(kw.*r).^2).^(11/6); % Von Karman Vertical Spectrum
        Sw = Sww;
        Sww=sqrt(Sww);
        Sww(1+TPointB/2:TPointB)=0;    
        dk = kw(1);

        ExpPart=zeros(1,TPointB);
        ExpPart(1:TPointB/2)=2.*pi.*rand(TPointB/2,1); 

        w=ifft(sqrt(2.*dk*r/Tw)*Sww.*exp(1j.*ExpPart));
        w=TPointB.*real(w.*exp(1j.*(1:1:TPointB).*2.*pi/TPointB));

        k1=2*pi*kw./2;       % Reduced velocity in Sears incl. 2pi fact and working with b=B/2
        Sears=0.5*0.13^2./(0.13^2+k1.^2)+0.5./(1+k1.^2)-1j.*k1.*(0.5*0.13./(0.13^2+k1.^2)+0.5./(1+k1.^2));
        Sears(1+TPointB/2:TPointB)=flip(conj(Sears));

        w_conv=real(ifft(fft(w).*Sears))'; %Incl. Sears admittance
end