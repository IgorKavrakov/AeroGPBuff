function Example1d_Plots(GP,FP)
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

nby=1;
FigWidth=18.48; nbx=2;spacey=1; spacex=2.18; leftmargin=1.2; rightmargin=0.5; topmargin=0.2; bottommargin=0.9;
FigDepth=nby*4+topmargin+bottommargin+(nby-1)*spacey;
FontAxis=7; FontLabel=9; FontLegend=7;
set(groot,'defaultLineLineWidth',0.6)

%% Plot 1  Flutter time history plot
f=figure('visible','on');clf(f);set(gcf, 'PaperUnits', 'centimeters');set(gcf, 'PaperSize', [FigWidth FigDepth]);set(gcf, 'PaperPositionMode', 'manual');set(gcf, 'PaperPosition', [0 0 FigWidth FigDepth]);
[ positions ] = subplot_pos(FigWidth,FigDepth,leftmargin,rightmargin,bottommargin,topmargin,nbx,nby,spacex,spacey);

FP(1).Pred.tau=(0:FP(1).Pred.Samp)*FP(1).Par.ds;
ax=axes('position',positions{1},'Layer','top');hold on;
h3=plot(FP(1).Pred.tau,FP(5).Pred.u(:,2)*180/pi,'-','color',Colors(3,:));
h2=plot(FP(1).Pred.tau,FP(4).Pred.u(:,2)*180/pi,'-','color',Colors(4,:));
h1=plot(FP(1).Pred.tau,FP(3).Pred.u(:,2)*180/pi,'-','color',Colors(1,:));
set(gca, 'FontSize',FontAxis);
xlabel('$\tau$ [-]','Interpreter', 'latex','FontSize',FontLabel)
ylabel('$\alpha_a$ [deg]','Interpreter', 'latex','FontSize',FontLabel)
ylim([-0.8 0.8]);yticks([-0.8:0.2:0.8]);ytickformat('%.1f');%xtickformat('%.1f');
l=legend([h1,h2,h3],{'Analytical Damped $U/(Bf_m)$=13.11','Analytical Critical $U_{cr}/(Bf_m)$=13.33','Analytical Divergent $U/(Bf_m)$=13.40'},'location','best');
set(l,'FontSize',FontLegend,'Interpreter', 'latex');
grid on; box on;

ax=axes('position',positions{2},'Layer','top');hold on;
h3=plot(FP(1).Pred.tau,GP(3).u(:,2)*180/pi,'-','color',Colors(3,:));
h2=plot(FP(1).Pred.tau,GP(2).u(:,2)*180/pi,'-','color',Colors(4,:));
h1=plot(FP(1).Pred.tau,GP(1).u(:,2)*180/pi,'-','color',Colors(1,:));
set(gca, 'FontSize',FontAxis);
xlabel('$\tau$ [-]','Interpreter', 'latex','FontSize',FontLabel)
ylabel('$\alpha_a$ [deg]','Interpreter', 'latex','FontSize',FontLabel)
ylim([-0.8 0.8]);yticks([-0.8:0.2:0.8]);ytickformat('%.1f');%xtickformat('%.1f');
l=legend([h1,h2,h3],{'GP Damped $U/(Bf_m)$=12.90','GP Critical $U_{cr}/(Bf_m)$=13.03','GP Divergent $U/(Bf_m)$=13.11'},'location','best');
set(l,'FontSize',FontLegend,'Interpreter', 'latex');
grid on; box on;

print(gcf,'Example1_FlatPlateAnalytical/Ex1d_Flut','-dpdf')


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