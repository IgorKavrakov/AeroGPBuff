%% by Igor Kavrakov (ik380@cam.ac.uk; igor.kavrakov@gmail.com)
clear all;  clc;  restoredefaultpath; matlabrc; close all;
addpath(genpath('GP'),'Example1_FlatPlateAnalytical');

% AeroGPBuff: Aeroelastic Analyses of Structures in Turbulent Wind Conditions using Gaussian Processes
% Please cite our work when you are you are using our software in your research or publications:

% Kavrakov, I., Morgenthal, G., and McRobie, A. 2024. Data-driven Aeroelastic Analyses of Structures in Turbulent Wind Conditions using Enhanced Gaussian Processes with Aerodynamic Priors. 
% J. Wind Eng. Ind. Aerodyn., 222, 104911. 
% https://doi.org/10.1016/j.jweia.2022.104911

% Accepted manuscript on arXiv:
% https://arxiv.org/abs/2406.15603

% The script includes a forced-vibration analysis with free-stream turbulence for the flat plate example for using a GP model.
% The script is based on the aforementioned article
% This script includes the results from Sec. 3 (Fundamental Application: Flat Plate)

%%%%%%%%% COPYRIGHT NOTICE %%%%%%%%% 
%  This file is part of AeroGPBuff.
%  AeroGPBuff is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation, either version 3 of the License, or
%  (at your option) any later version.
% 
%  AeroGPBuff is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
% 
%  You should have received a copy of the GNU General Public License
%  along with AeroGPBuff.  If not, see <https://www.gnu.org/licenses/>.

% Copyright (c) Igor Kavrakov, Guido Morgenthal, Allan McRobie 2024
fprintf(['AeroGPBuff: Aeroelastic Analyses of Structures in Turbulent Wind Conditions using Gaussian Processes \nIgor Kavrakov, Guido Morgenthal, Allan McRobie 2024 (c) \nCite as:\n Kavrakov, I., McRobie, A., and Morgenthal, G. 2022.\n Data-driven aerodynamic analysis of structures using Gaussian Processes.\n J. Wind Eng. Ind. Aerodyn., 222, 104911.\n https://doi.org/10.1016/j.jweia.2022.104911\n\n']);

%% Control
FP.Par.ds=0.05;            %Reduced time-step (Warning: be consistent with the training inpu from Example 1a_FlatPlateForced.m)

%Dynamic properties for flutter analysis
FP.Pred.Excitation='Aeroelastic'; %Aeroelastic simulation - Flutter
FP.Par.B=31;               % Width [m]
FP.Par.m=22.7400*1000;  % Mass per unit length [kg/m]
FP.Par.I=2470*1000;     % Mass moment of inertia [kgm^2/m]
FP.Par.fh=0.1;          % Vertical frequency [Hz]
FP.Par.fa=0.278;          % Torsional frequency [Hz]
FP.Par.psi=0.003;       % Damping ratio [-]
FP.Par.Redtime=400;     % Reduced time for analysis [-]
FP.Par.rho=1.2;         % Density [kg/m^3]
FP.Par.U_r=[75.58,...   % Test flutter velocities [m/s]
            76.3,...    % Reduced velocity: U_r/(B*(fh+fa)/2): 
            76.8,...    % Reduced velocity: U_r/(B*(fh+fa)/2):                        
            78.1005,... % [12.90 13.03 13.11 13.33 13.40] [-]
            78.8035];  
        
%% Training - Based on Example1c
if isfile('Example1_FlatPlateAnalytical/Example1_Train.mat')
load('Example1_FlatPlateAnalytical/Example1_Train','GP');%Load hyperparameters 
else
error('Please supply the hyperparameters (see Example1c_FlatPlateAerodynamicPrior_Forced');    
end
%%  Aeroelastic prediction
%Get alpha - it is not needed to be inverted every step.
    [~,~,~,~,~,~,~,GP.Train(1).alpha]=GP_Process(GP.Train(1).x,GP.Train(1).y,GP.Train(1).x(1,:),GP.Train(1).hyp,GP.Par.jitter,GP.Par.Noise,GP.meanFlag(1)); 
    [~,~,~,~,~,~,~,GP.Train(2).alpha]=GP_Process(GP.Train(2).x,GP.Train(2).y,GP.Train(2).x(1,:),GP.Train(2).hyp,GP.Par.jitter,GP.Par.Noise,GP.meanFlag(2)); 
%% Flutter loop
for j=1:length(FP(1).Par.U_r)
    FP(j)=FP(1); 
    FP(j).Par.U=FP(j).Par.U_r(j);
    FP(j)=Example1_FlatPlateAnalytical(FP(j),GP(1).Par.Lag); %FP Analytical prediction 

    GP(j).meanFlag=GP(1).meanFlag;
    GP(j).Train=GP(1).Train;
    GP(j).Pred(1).y=zeros(FP(1).Pred.Samp,1); %CL
    GP(j).Pred(2).y=zeros(FP(1).Pred.Samp,1); %CM
    GP(j).u=zeros(FP(1).Pred.Samp,2); %Dispalcements
    GP(j).u_d=zeros(FP(1).Pred.Samp,2); %Velocity
    GP(j).u_2d=zeros(FP(1).Pred.Samp,2); %Acceleration
    GP(j).u(1,1)=0.5; %Initial conditions

    alpha_h=zeros(1,1+GP(1).Par.Lag);  alpha_a=zeros(1,1+GP(1).Par.Lag);  alpha_h_d=0; alpha_a_d=0;
    for i=1:FP(j).Pred.Samp
        
        alpha_h(1,1)=atan(GP(j).u_d(i,1)/FP(1).Par.B);  
        alpha_a(1,1)=GP(j).u(i,2);              
        alpha_h_d   =1/(1+tan(alpha_h(1,1))^2)*(GP(j).u_2d(i,1)/FP(1).Par.B); 
        alpha_a_d   =GP(j).u_d(i,2); 

        %Calculate forces
        GP(j).Pred(1).x=[alpha_h alpha_a alpha_h_d alpha_a_d alpha_h*0 alpha_h_d*0]; 
        GP(j).Pred(1).x=GP(j).Pred(1).x./GP(j).Train(1).x_max;
        GP(j).Pred(1).y(i)=GP_Process_Predict_Mean(GP(j).Train(1).x,GP(j).Train(1).y,GP(j).Pred(1).x,GP(j).Train(1).hyp,GP(j).Train(1).alpha,GP(j).meanFlag(1)); %Predictions

        GP(j).Pred(2).x=[alpha_h alpha_a alpha_h_d alpha_a_d  alpha_h*0 alpha_h_d*0]; 
        GP(j).Pred(2).x=GP(j).Pred(2).x./GP(j).Train(2).x_max;
        GP(j).Pred(2).y(i)=GP_Process_Predict_Mean(GP(j).Train(2).x,GP(j).Train(2).y,GP(j).Pred(2).x,GP(j).Train(2).hyp,GP(j).Train(2).alpha,GP(j).meanFlag(2)); %Predictions
        
        p=FP(j).Par.Pres.*[GP(j).Pred(1).y(i).*GP(j).Train(1).y_max GP(j).Pred(2).y(i).*GP(j).Train(2).y_max];
        [GP(j).u(i+1,:),GP(j).u_d(i+1,:),GP(j).u_2d(i+1,:)] = NewmarkSDOF(p,GP(j).u(i,:),GP(j).u_d(i,:),GP(j).u_2d(i,:),1/4,1/2,FP(1).Par.ds,FP(j).Par.M,FP(j).Par.C,FP(j).Par.K);    

        for ll=1:GP(1).Par.Lag % Add the lags
            alpha_h(1,end+1-ll)=alpha_h(1,end-ll);
            alpha_a(1,end+1-ll)=alpha_a(1,end-ll);              
        end

    end

end

save('Example1_FlatPlateAnalytical/Ex1d_Out','GP','FP','-v7.3');
%% Plots
load('Example1_FlatPlateAnalytical/Ex1d_Out','GP','FP');
Example1d_Plots(GP,FP)
