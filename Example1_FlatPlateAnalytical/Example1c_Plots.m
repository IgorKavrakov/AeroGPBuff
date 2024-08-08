function Example1c_Plots(GP,FP)
% This is a speciic function for plots in the Flat Plate Example

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

%% Plot properties
close all;
Colors=[57  106 177
        218 124 48
        62  150 81
        204 37  41
        83  81  84
        107 76  154
        146 36  40
        148 139 61]./255;


FigWidth=18.48;nby=2; nbx=3; spacey=1; spacex=1.0; leftmargin=1.2; rightmargin=0.5; topmargin=0.2; bottommargin=0.9;
FigDepth=nby*4+topmargin+bottommargin+(nby-1)*spacey;
FontAxis=7; FontLabel=9; FontLegend=7;
set(groot,'defaultLineLineWidth',0.6)

%% Plot 1 - Training set input
TW = [400 500];
f=figure('visible','on');clf(f);set(gcf, 'PaperUnits', 'centimeters');set(gcf, 'PaperSize', [FigWidth FigDepth]);set(gcf, 'PaperPositionMode', 'manual');set(gcf, 'PaperPosition', [0 0 FigWidth FigDepth]);
[ positions ] = subplot_pos(FigWidth,FigDepth,leftmargin,rightmargin,bottommargin,topmargin,nbx,nby,spacex,spacey);

GP.Train(1).x=GP.Train(1).x*GP.Train(1).x_max;% Denormalize
GP.Train(2).x=GP.Train(2).x*GP.Train(2).x_max;% Denormalize

ax=axes('position',positions{4},'Layer','top');hold on;
h10=fill([TW(1),TW(2),TW(2),TW(1)]',[-10,-10,10,10]',[0.8 0.8 0.8],'LineStyle','none','facealpha',.5);
plot(FP.Train.tau,GP.Train(1).x(:,1)*180/pi,'-','color',Colors(1,:));
set(gca, 'FontSize',FontAxis);
xlabel('$\tau$ [-]','Interpreter', 'latex','FontSize',FontLabel)
ylabel('$\alpha_h$ [deg]','Interpreter', 'latex','FontSize',FontLabel)
% ytickformat('%.1f');
% xlim([0 280]);xticks(0:40:280);
grid on; box on;

ax=axes('position',positions{5},'Layer','top');hold on;
h10=fill([TW(1),TW(2),TW(2),TW(1)]',[-10,-10,10,10]',[0.8 0.8 0.8],'LineStyle','none','facealpha',.5);
plot(FP.Train.tau,GP.Train(1).x(:,2+GP.Par.Lag)*180/pi,'-','color',Colors(4,:));
set(gca, 'FontSize',FontAxis);
xlabel('$\tau$ [-]','Interpreter', 'latex','FontSize',FontLabel)
ylabel('$\alpha_a$ [deg]','Interpreter', 'latex','FontSize',FontLabel)
grid on; box on;

ax=axes('position',positions{6},'Layer','top');hold on;
h10=fill([TW(1),TW(2),TW(2),TW(1)]',[-15,-15,15,15]',[0.8 0.8 0.8],'LineStyle','none','facealpha',.5);
plot(FP.Train.tau,GP.Train(1).x(:,5+2*GP.Par.Lag)*180/pi,'-','color',Colors(3,:));
set(gca, 'FontSize',FontAxis);
xlabel('$\tau$ [-]','Interpreter', 'latex','FontSize',FontLabel)
ylabel('$\alpha_w$ [deg]','Interpreter', 'latex','FontSize',FontLabel)
grid on; box on;



FigWidth=18.48;nby=2; nbx=2;spacey=1; spacex=2.18; leftmargin=1.2; rightmargin=0.5; topmargin=0.2; bottommargin=0.9;
FigDepth=nby*4+topmargin+bottommargin+(nby-1)*spacey;
FontAxis=7; FontLabel=9; FontLegend=7;
[ positions ] = subplot_pos(FigWidth,FigDepth,leftmargin,rightmargin,bottommargin,topmargin,nbx,nby,spacex,spacey);

GP.Pred(1).y=GP.Pred(1).y*GP.Train(1).y_max;% Denormalize
GP.Pred(2).y=GP.Pred(2).y*GP.Train(2).y_max;% Denormalize
GP.Train(1).y=GP.Train(1).y*GP.Train(1).y_max;% Denormalize
GP.Train(2).y=GP.Train(2).y*GP.Train(2).y_max;% Denormalize
GP.Pred(1).y_targ=GP.Pred(1).y_targ*GP.Train(1).y_max;% Denormalize
GP.Pred(2).y_targ=GP.Pred(2).y_targ*GP.Train(2).y_max;% Denormalize
GP.Pred(1).y_std=sqrt(diag(GP.Pred(1).y_var))*GP.Train(1).y_max; %Denormalize
GP.Pred(2).y_std=sqrt(diag(GP.Pred(2).y_var))*GP.Train(2).y_max; %Denormalize
GP.Pred(1).m=GP.Pred(1).m*GP.Train(1).y_max;% Denormalize
GP.Pred(2).m=GP.Pred(2).m*GP.Train(2).y_max;% Denormalize

ax=axes('position',positions{1},'Layer','top');hold on;
% h10=fill([FP.Train.tau; flip(FP.Train.tau)], [GP.Pred(1).y+2.56.*GP.Pred(1).y_std; flip(GP.Pred(1).y-2.56*GP.Pred(1).y_std)],Colors(1,:),'LineStyle','none','facealpha',.3);
h4=plot(FP.Train.tau,GP.Pred(1).m,'--','color',Colors(3,:));
h1=plot(FP.Train.tau,GP.Pred(1).y_targ,'-','color',Colors(4,:));
h3=plot(FP.Train.tau,GP.Pred(1).y,'--','color',Colors(1,:),'LineWidth',0.8);
set(gca, 'FontSize',FontAxis);
xlabel('$\tau$ [-]','Interpreter', 'latex','FontSize',FontLabel)
ylabel('$C_L$ [-]','Interpreter', 'latex','FontSize',FontLabel)
% ytickformat('%.2f');%xtickformat('%.1f');
l=legend([h1 h3 h4],{'Analytical','GP Train','QS Prior'},'location','best');
set(l,'FontSize',FontLegend,'Interpreter', 'latex');
xlim([TW(1) TW(2)])
grid on; box on;

ax=axes('position',positions{2},'Layer','top');hold on;
h10=fill([FP.Train.tau; flip(FP.Train.tau)], [GP.Pred(2).y+2.56.*GP.Pred(2).y_std; flip(GP.Pred(2).y-2.56*GP.Pred(2).y_std)],Colors(1,:),'LineStyle','none','facealpha',.3);
h4=plot(FP.Train.tau,GP.Pred(2).m,'--','color',Colors(3,:));
h1=plot(FP.Train.tau,GP.Pred(2).y_targ,'-','color',Colors(4,:));
h3=plot(FP.Train.tau,GP.Pred(2).y,'--','color',Colors(1,:),'LineWidth',0.8);
set(gca, 'FontSize',FontAxis);
xlabel('$\tau$ [-]','Interpreter', 'latex','FontSize',FontLabel)
ylabel('$C_M$ [-]','Interpreter', 'latex','FontSize',FontLabel)
ytickformat('%.1f');%xtickformat('%.1f');
l=legend([h1 h3 h4],{'Analytical','GP Train','QS Prior'},'location','best');
set(l,'FontSize',FontLegend,'Interpreter', 'latex');
xlim([TW(1) TW(2)])
grid on; box on;
f.Renderer='Painters';
print(gcf,'Example1_FlatPlateAnalytical/Ex1c_Train','-dpdf')

%% Plot 2 - Flutter Derivatives
%Preprocess Lift and moment due to alpha_h
GP.Pred(4).x=GP.Pred(4).x*GP.Train(4).x_max;% Denormalize
GP.Pred(4).y=GP.Pred(4).y*GP.Train(4).y_max;% Denormalize
GP.Pred(4).y_targ=GP.Pred(4).y_targ*GP.Train(4).y_max;% Denormalize

GP.Pred(5).x=GP.Pred(5).x*GP.Train(5).x_max;% Denormalize
GP.Pred(5).y=GP.Pred(5).y*GP.Train(5).y_max;% Denormalize
GP.Pred(5).y_targ=GP.Pred(5).y_targ*GP.Train(5).y_max;% Denormalize

GP.Pred(6).x=GP.Pred(6).x*GP.Train(6).x_max;% Denormalize
GP.Pred(6).y=GP.Pred(6).y*GP.Train(6).y_max;% Denormalize
GP.Pred(6).y_targ=GP.Pred(6).y_targ*GP.Train(6).y_max;% Denormalize

GP.Pred(7).x=GP.Pred(7).x*GP.Train(5).x_max;% Denormalize
GP.Pred(7).y=GP.Pred(7).y*GP.Train(5).y_max;% Denormalize
GP.Pred(7).y_targ=GP.Pred(7).y_targ*GP.Train(7).y_max;% Denormalize

StepsVr=cumsum([1 FP.Pred.Cycl.*FP.Pred.Vr./FP.Par.ds+1]);

for i=1:length(FP.Pred.Vr)
 S=StepsVr(i);
 E=StepsVr(i+1)-1;
 
    [H1(i),H4(i)] = DerFit(FP.Pred.Vr(i),atan(GP.Pred(4).x(S:E,1)),atan(GP.Pred(4).x(S:E,3+2*GP.Par.Lag)),GP.Pred(4).y(S:E),1);
    [A1(i),A4(i)] = DerFit(FP.Pred.Vr(i),atan(GP.Pred(5).x(S:E,1)),atan(GP.Pred(5).x(S:E,3+2*GP.Par.Lag)),GP.Pred(5).y(S:E),1);
    [H2(i),H3(i)] = DerFit(FP.Pred.Vr(i),GP.Pred(6).x(S:E,2+GP.Par.Lag),GP.Pred(6).x(S:E,4+2*GP.Par.Lag),GP.Pred(6).y(S:E),0);
    [A2(i),A3(i)] = DerFit(FP.Pred.Vr(i),GP.Pred(7).x(S:E,2+GP.Par.Lag),GP.Pred(7).x(S:E,4+2*GP.Par.Lag),GP.Pred(7).y(S:E),0);
end
K=2*pi./FP.Pred.Vr;
C=conj(besselh(1,K/2))./(conj(besselh(1,K/2))+ 1j.*conj(besselh(0,K/2)));
G=imag(C); F=real(C);

H1_ana=-2.*pi.*F./K;              
H2_ana=-pi./2*(1+4.*G./K+F)./K;     
H3_ana=-pi.*(2.*F-0.5.*G.*K)./(K.*K); 
H4_ana=pi.*0.5.*(1+4.*G./K);        
A1_ana=0.5.*pi.*F./K;              
A2_ana=-pi./(2.*K.*K).*(K./4-G-K.*F./4); 
A3_ana=0.5.*pi.*(K.*K./32+F-K.*G./4)./(K.*K); 
A4_ana=-0.5.*pi.*G./K; 

% Plot 2 subplots
nby=1;
FigDepth=nby*4+topmargin+bottommargin+(nby-1)*spacey;
f=figure('visible','on');clf(f);set(gcf, 'PaperUnits', 'centimeters');set(gcf, 'PaperSize', [FigWidth FigDepth]);set(gcf, 'PaperPositionMode', 'manual');set(gcf, 'PaperPosition', [0 0 FigWidth FigDepth]);
[ positions ] = subplot_pos(FigWidth,FigDepth,leftmargin,rightmargin,bottommargin,topmargin,nbx,nby,spacex,spacey);

ax=axes('position',positions{1,1},'Layer','top');hold on;
h1=plot(FP.Pred.Vr,H1_ana,'-','color',Colors(1,:));
plot(FP.Pred.Vr,H1,'.-','color',Colors(1,:));
h2=plot(FP.Pred.Vr,H2_ana,'-','color','k');
plot(FP.Pred.Vr,H2,'.-','color','k');
h3=plot(FP.Pred.Vr,H3_ana/2,'-','color',Colors(4,:));
plot(FP.Pred.Vr,H3/2,'.-','color',Colors(4,:));
h4=plot(FP.Pred.Vr,H4_ana,'-','color',Colors(3,:));
plot(FP.Pred.Vr,H4,'.-','color',Colors(3,:));

set(gca, 'FontSize',FontAxis);
xlabel('$V_r$ [-]','Interpreter', 'latex','FontSize',FontLabel)
ylabel('$H_1^*,H_2^*,H_3^*/2,H_4^*$ [-]','Interpreter', 'latex','FontSize',FontLabel)
ylimit=get(gca,'ylim'); plot([14 14],ylimit,'k'); plot([2 2],ylimit,'k');
l=legend([h1 h2 h3 h4],{'$H_1^*$','$H_2^*$','$H_3^*$','$H_4^*$'},'location','southwest');
set(l,'FontSize',FontLegend,'Interpreter', 'latex');
xlim([0 16]);xticks(0:2:16);
grid on; box on;

ax=axes('position',positions{2,1},'Layer','top');hold on;
h1=plot(FP.Pred.Vr,A1_ana,'-','color',Colors(1,:));
plot(FP.Pred.Vr,A1,'.-','color',Colors(1,:));
h2=plot(FP.Pred.Vr,A2_ana,'-','color','k');
plot(FP.Pred.Vr,A2,'.-','color','k');
h3=plot(FP.Pred.Vr,A3_ana/2,'-','color',Colors(4,:));
plot(FP.Pred.Vr,A3/2,'.-','color',Colors(4,:));
h4=plot(FP.Pred.Vr,A4_ana,'-','color',Colors(3,:));
plot(FP.Pred.Vr,A4,'.-','color',Colors(3,:));

set(gca, 'FontSize',FontAxis);
xlabel('$V_r$ [-]','Interpreter', 'latex','FontSize',FontLabel)
ylabel('$A_1^*,A_2^*,A_3^*/2,A_4^*$ [-]','Interpreter', 'latex','FontSize',FontLabel)
ylim([-4 6])
ylimit=get(gca,'ylim'); plot([14 14],ylimit,'k'); plot([2 2],ylimit,'k');
l=legend([h1 h2 h3 h4],{'$A_1^*$','$A_2^*$','$A_3^*$','$A_4^*$'},'location','northwest');
set(l,'FontSize',FontLegend,'Interpreter', 'latex');
xlim([0 16]);xticks(0:2:16);
grid on; box on;
print(gcf,'Example1_FlatPlateAnalytical/Ex1c_FlutterDer','-dpdf')

%% Plot 3 - Static wind coefficients & Admittance
%Preprocess Lift and moment due to alpha
GP.Pred(8).x=GP.Pred(8).x*GP.Train(8).x_max;% Denormalize
GP.Pred(8).y=GP.Pred(8).y*GP.Train(8).y_max;% Denormalize
GP.Pred(8).y_targ=GP.Pred(8).y_targ*GP.Train(8).y_max;% Denormalize

GP.Pred(9).x=GP.Pred(9).x*GP.Train(9).x_max;% Denormalize
GP.Pred(9).y=GP.Pred(9).y*GP.Train(9).y_max;% Denormalize
GP.Pred(9).y_targ=GP.Pred(9).y_targ*GP.Train(9).y_max;% Denormalize

FP.Pred.StatTau=20;          %P Cycles - number of cycles for prediction
Inc=(FP.Pred.StatTau./FP.Par.ds+1); %Initiate

CL=zeros(length(FP.Pred.StatAlpha),1);
CM=zeros(length(FP.Pred.StatAlpha),1);
for i=1:length(FP.Pred.StatAlpha)
CL(i)=mean(GP.Pred(8).y(1+(i-1)*Inc:i*Inc))*-1;  %-4 as CM=CL/-4
CM(i)=mean(GP.Pred(9).y(1+(i-1)*Inc:i*Inc))*4;     
end
C_Analytical=FP.Pred.StatAlpha*pi/180*pi*2; %Analytical CL

%Preprocess Lift and moment due to alpha_h
GP.Pred(10).x=GP.Pred(10).x*GP.Train(10).x_max;% Denormalize
GP.Pred(10).y=GP.Pred(10).y*GP.Train(10).y_max;% Denormalize
GP.Pred(10).y_targ=GP.Pred(10).y_targ*GP.Train(10).y_max;% Denormalize

GP.Pred(11).x=GP.Pred(11).x*GP.Train(11).x_max;% Denormalize
GP.Pred(11).y=GP.Pred(11).y*GP.Train(11).y_max;% Denormalize
GP.Pred(11).y_targ=GP.Pred(11).y_targ*GP.Train(11).y_max;% Denormalize
StepsVr=cumsum([1 FP.Pred.CyclW.*FP.Pred.VrW./FP.Par.ds+1]);

for i=1:length(FP.Pred.VrW)
 S=StepsVr(i);
 E=StepsVr(i+1)-1;
 
    w_U0=FFT(tan(GP.Pred(10).x(S:E,5+2*GP.Par.Lag)));

    CL0 =FFT(GP.Pred(10).y(S:E));
    phiCL=PhaseDifference(GP.Pred(10).x(S:E,5+2*GP.Par.Lag),GP.Pred(10).y(S:E));
    Signal=CL0*exp(1j*phiCL)./(-2*pi*w_U0); %check signs for 1j
    ImagCL(i)=imag(Signal);
    RealCL(i)=real(Signal);
    
    CM0 =FFT(GP.Pred(11).y(S:E));
    phiCM=PhaseDifference(GP.Pred(11).x(S:E,5+2*GP.Par.Lag),GP.Pred(11).y(S:E));
    Signal=CM0*exp(1j*phiCM)./(pi/2*w_U0); %check signs for 1j
    ImagCM(i)=imag(Signal);
    RealCM(i)=real(Signal);
    
end
k1=pi./FP.Pred.VrW;
Sears=0.5*0.13^2./(0.13^2+k1.^2)+0.5./(1+k1.^2)-1j.*k1.*(0.5*0.13./(0.13^2+k1.^2)+0.5./(1+k1.^2));

% Plot 2 subplots
nby=1;
FigDepth=nby*4+topmargin+bottommargin+(nby-1)*spacey;
f=figure('visible','on');clf(f);set(gcf, 'PaperUnits', 'centimeters');set(gcf, 'PaperSize', [FigWidth FigDepth]);set(gcf, 'PaperPositionMode', 'manual');set(gcf, 'PaperPosition', [0 0 FigWidth FigDepth]);
[ positions ] = subplot_pos(FigWidth,FigDepth,leftmargin,rightmargin,bottommargin,topmargin,nbx,nby,spacex,spacey);

ax=axes('position',positions{1},'Layer','top');hold on;
h1=plot(FP.Pred.StatAlpha,CL,'.-','color',Colors(1,:));
h2=plot(FP.Pred.StatAlpha,CM,'.-','color',Colors(4,:));
h3=plot(FP.Pred.StatAlpha,C_Analytical,'-','color','k');
set(gca, 'FontSize',FontAxis);
xlabel('$\alpha$ [deg]','Interpreter', 'latex','FontSize',FontLabel)
ylabel('$-\tilde{C}_L,4\tilde{C}_M$ [-]','Interpreter', 'latex','FontSize',FontLabel)
l=legend([h1 h2 h3],{'GP $\tilde{C}_L$','GP $\tilde{C}_M$','Analytical'},'location','northwest');
set(l,'FontSize',FontLegend,'Interpreter', 'latex');
xlim([-1 1]);xticks(-1:0.5:1);
ylim([-0.15 0.15]);yticks(-0.15:0.05:0.15);
xtickformat('%.1f');ytickformat('%.2f');
grid on; box on;

ax=axes('position',positions{2},'Layer','top');hold on;
h1=plot(FP.Pred.VrW,real(Sears),'-','color','k');
h2=plot(FP.Pred.VrW,RealCL,'.-','color',Colors(1,:));
h3=plot(FP.Pred.VrW,RealCM,'.-','color',Colors(4,:));
h4=plot(FP.Pred.VrW,imag(Sears),'--','color','k');
h5=plot(FP.Pred.VrW,ImagCL,'.--','color',Colors(1,:));
h6=plot(FP.Pred.VrW,ImagCM,'.--','color',Colors(4,:));
set(gca, 'FontSize',FontAxis);
xlabel('$V_r$ [-]','Interpreter', 'latex','FontSize',FontLabel)
ylabel('$\chi_L,\chi_M$ [-]','Interpreter', 'latex','FontSize',FontLabel)
ylim([-0.4 1.0]);yticks([-0.4:0.2:1.0]);ytickformat('%.1f');%xtickformat('%.1f');
ylimit=get(gca,'ylim'); plot([2 2],ylimit,'k');%  plot([14 14],ylimit,'k');
l=legend([h1 h2 h3 ],{'Analytical - Sears','GP $\chi_L$','GP $\chi_M$'},'location','east');
set(l,'FontSize',FontLegend,'Interpreter', 'latex');
xlim([0 30]);xticks(0:5:30);
grid on; box on;

print(gcf,'Example1_FlatPlateAnalytical/Ex1c_SWCAdm','-dpdf')

end
%% Help functions
%Plot function
function [ positions ] = subplot_pos(plotwidth,plotheight,leftmargin,rightmargin,bottommargin,topmargin,nbx,nby,spacex,spacey)
    subxsize=(plotwidth-leftmargin-rightmargin-spacex*(nbx-1.0))/nbx;
    subysize=(plotheight-topmargin-bottommargin-spacey*(nby-1.0))/nby;
    for i=1:nbx
       for j=1:nby
           xfirst=leftmargin+(i-1.0)*(subxsize+spacex);
           yfirst=bottommargin+(j-1.0)*(subysize+spacey);
           positions{i,j}=[xfirst/plotwidth yfirst/plotheight subxsize/plotwidth subysize/plotheight];
       end
    end
end

%Derivative fit function
function [Der1,Der2] = DerFit(Vr,alpha,alpha_d,F,Type)
K=2*pi/Vr;
if Type %Vertical Displacement Related Derivatives
    Ss(:,1)=K.*alpha;
    Ss(:,2)=-alpha_d;
else    %Rotation Related Derivatives
    Ss(:,1)=K*alpha_d;
    Ss(:,2)=K.^2.*alpha;
end
StF=Ss'*F; %Perform least square fit
StS=Ss'*Ss;
Der=StS\StF;

Der1=Der(1); %Velocity related derivative
Der2=Der(2); %Displacement related derivative
end

%Admittance fit function
function [X_max] = FFT(signal)
%% By Igor Kavrakov
signal(:)=signal;
if mod(length(signal),2)==1
    signal=signal(1:end-1);
end

NSign=length(signal);
X=fft(signal); %X[k]=sum_n_Nsignal(x[n]*W_Nsignal^{kn})=sum_n_Nsignal(x[n]*(-j)^nk) %TWO SIDED!
X=X./NSign; 

X=X(1:ceil((NSign+1)/2));

X=abs(X);
X(1:end-1)=X(1:end-1)*2;
X_max = max(X);
end

function phi = PhaseDifference(x,y)
% Phase Difference Measurement  
if size(x,2)>1
    x = x';
end
% represent y as column-vector if it is not
if size(y,2)>1
    y = y';
end
% remove the DC component
x = x-mean(x);
y = y-mean(y);

X = fft(x);

Y = fft(y);

[~,indx] = max(abs(X));
[~,indy] = max(abs(Y));

px = angle(X(indx));
py = angle(Y(indy));

phi=py-px;
end