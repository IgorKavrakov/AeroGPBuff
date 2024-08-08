function Example1e_Plots(GP,FP)
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


FigWidth=18.48;nby=2; nbx=2;spacey=1; spacex=2.18; leftmargin=1.2; rightmargin=0.5; topmargin=0.2; bottommargin=0.9;
FigDepth=nby*4+topmargin+bottommargin+(nby-1)*spacey;
FontAxis=7; FontLabel=9; FontLegend=7;
set(groot,'defaultLineLineWidth',0.6)

%% Plot 1 STD of the buffeting response
f=figure('visible','on');clf(f);set(gcf, 'PaperUnits', 'centimeters');set(gcf, 'PaperSize', [FigWidth FigDepth]);set(gcf, 'PaperPositionMode', 'manual');set(gcf, 'PaperPosition', [0 0 FigWidth FigDepth]);
[ positions ] = subplot_pos(FigWidth,FigDepth,leftmargin,rightmargin,bottommargin,topmargin,nbx,nby,spacex,spacey);
  
% Preprocess to find the STD (i.e. RMS without the mean)
   U_red = FP(1).Par.U_r / FP(1).Par.B / (FP(1).Par.fh*0.5  + FP(1).Par.fa*0.5);
for j=1:length(FP(1).Par.U_r)
   h_rms(j,1)=std(atan(FP(j).Pred.u_d(:,1)./FP(1).Par.B)*180/pi);
   h_rms(j,2)=std(atan(GP(j).u_d(:,1)/FP(1).Par.B)*180/pi);   
   a_rms(j,1)=std(FP(j).Pred.u(:,2)*180/pi);
   a_rms(j,2)=std(GP(j).u(:,2)*180/pi);   
end



ax=axes('position',positions{1,1},'Layer','top');hold on;
h2=plot(U_red,h_rms(:,2),'o-','markersize',3,'color',Colors(1,:));
h1=plot(U_red,h_rms(:,1),'.-','color',Colors(4,:));
set(gca, 'FontSize',FontAxis);
xlabel('$U/(Bf_m)$ [-]','Interpreter', 'latex','FontSize',FontLabel)
ylabel('$\sigma_{\alpha_h}$ [deg]','Interpreter', 'latex','FontSize',FontLabel)
xlim([3.0 12.0]);xticks([3.0:1:12.0]);%xtickformat('%.1f');%xtickformat('%.1f');
ylim([0 2]);yticks([0:0.25:2]);ytickformat('%.2f');%xtickformat('%.1f');
l=legend([h1,h2],{'Analytical','GP'},'location','best');
set(l,'FontSize',FontLegend,'Interpreter', 'latex');
grid on; box on;

ax=axes('position',positions{2,1},'Layer','top');hold on;
h2=plot(U_red,a_rms(:,2),'o-','markersize',3,'color',Colors(1,:));
h1=plot(U_red,a_rms(:,1),'.-','color',Colors(4,:));
set(gca, 'FontSize',FontAxis);
xlabel('$U/(Bf_m)$ [-]','Interpreter', 'latex','FontSize',FontLabel)
ylabel('$\sigma_{\alpha_a}$ [deg]','Interpreter', 'latex','FontSize',FontLabel)
xlim([3 12]);xticks([3:1:12]);%xtickformat('%.1f');%xtickformat('%.1f');
ylim([0 4]);yticks([0.0:0.5:4]);ytickformat('%.1f');%xtickformat('%.1f');
set(l,'FontSize',FontLegend,'Interpreter', 'latex');
grid on; box on;

%% Plot 2  Time History at selected wind speed (j=5)
j = 5;
FP(1).Pred.tau=(0:FP(1).Pred.Samp)*FP(1).Par.ds;

ax=axes('position',positions{1,2},'Layer','top');hold on;
h1=plot(FP(1).Pred.tau,atan(FP(j).Pred.u_d(:,1)/FP(1).Par.B)*180/pi,'-','color',Colors(4,:));
h2=plot(FP(1).Pred.tau,atan(GP(j).u_d(:,1)/FP(1).Par.B)*180/pi,'--','color',Colors(1,:),'LineWidth',0.8);
set(gca, 'FontSize',FontAxis);
xlabel('$\tau$ [-]','Interpreter', 'latex','FontSize',FontLabel)
ylabel('$\alpha_h$ [deg]','Interpreter', 'latex','FontSize',FontLabel)
xlim([200 500])
l=legend([h1,h2],{['Analytical ($U/(Bf_m)$ = ' num2str(U_red(j),'%.1f') ')'],...
                  ['GP ($U/(Bf_m)$ = ' num2str(U_red(j),'%.1f') ')']},...
                  'location','best');
set(l,'FontSize',FontLegend,'Interpreter', 'latex');
grid on; box on;

ax=axes('position',positions{2,2},'Layer','top');hold on;
h1=plot(FP(1).Pred.tau,FP(j).Pred.u(:,2)*180/pi,'-','color',Colors(4,:));
h2=plot(FP(1).Pred.tau,GP(j).u(:,2)*180/pi,'--','color',Colors(1,:),'LineWidth',0.8);
set(gca, 'FontSize',FontAxis);
xlabel('$\tau$ [-]','Interpreter', 'latex','FontSize',FontLabel)
ylabel('$\alpha_a$ [deg]','Interpreter', 'latex','FontSize',FontLabel)
xlim([200 500])
set(l,'FontSize',FontLegend,'Interpreter', 'latex');
grid on; box on;

print(gcf,'Example1_FlatPlateAnalytical/Ex1e_Buf','-dpdf')


end

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