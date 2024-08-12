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
% 
% Copyright (c) Igor Kavrakov, Guido Morgenthal, Allan McRobie 2024
fprintf(['AeroGPBuff: Aeroelastic Analyses of Structures in Turbulent Wind Conditions using Gaussian Processes \nIgor Kavrakov, Guido Morgenthal, Allan McRobie 2024 (c) \nCite as:\n Kavrakov, I., McRobie, A., and Morgenthal, G. 2022.\n Data-driven aerodynamic analysis of structures using Gaussian Processes.\n J. Wind Eng. Ind. Aerodyn., 222, 104911.\n https://doi.org/10.1016/j.jweia.2022.104911\n\n']);

%% Control
%GP properties
GP.Par.Train=0;     %If Train=0 - no training; Hyperparameters are loaded (default as in Github)
GP.Par.Mean=1;      %Include mean based on linear quasi steady theory
GP.Par.Lag=200;     %Regressors for S_alpha - alpha,. Note if S_alpha=1, it means we take 1 previous i.e. Alpha(i-1).
GP.Par.jitter=eps;  %Jitter term
GP.Par.SNR=20;      %Signal to nose ratio
GP.Par.SubSet=3;    %Subset of regressors (F)
GP.Par.Noise=-1;    %Optimise for noise? -1 >yes
GP.Par.OptiControl1=optimset('GradObj','on',...     %Optimiser properties
                             'TolX',1e-8,...
                             'Display','iter',...
                             'MaxFunEval',200000,...
                             'MaxIter', 500);

%Flat plate input & training properties
FP.Par.B=31;                    % Cord Width [m]
FP.Train.Alpha_STD=[.5 2.5];    %Standard deviation (vector) of the effecitve angle - training
FP.Train.Alpha_AmpRel=.05;      %Relative Fourier amplitude - DoE for training (r_l)
FP.Train.Vr=[2 14];             %Input angle magnitude. Vr-Reduced velocity
FP.Train.Cycl=20;               %Cycles - number of cycles for training
FP.Par.ds=0.05;                 %Reduced time-step (Warning: use 0.05 i.e. be consistent if using FP_RandInp)

FP.Train.Buf=1;         %Training to include buffeting forces
FP.Par.Iw=[0.03 0.06];  %Turbulence intensity (vector)
FP.Par.r=2;             %Ratio lenght scale vs width Lw/B
FP.Par.Tw=2;            %Time scale of correlation (L.U) (e.g. 3 means it takes 3 seconds for a point to travel distance equivalent to average gust length scale Lw)
FP.Train.Vr_b=16;       %(Training) Gust largest reduced frequency to be the wind spectrum (this controls the length of time series) 
FP.Train.Cycl_b=10;     %(Training) Cycles
FP.Train.Coupled=1;     %Consider the same number of samples for buffeting & self-excited forces (based on the self-exccited forces)

%% Training
rng(1); %Reproducibility
if GP.Par.Train
    input('!!!!!!!!!!!!!!!!!!!!!!!!!! WARNING !!!!!!!!!!!!!!!!!!!!!!!!!!\n\nThe training of the GP model is computationally expensive (~16h - Linux server).\nFor prediction, use the determined hyperparameters by setting GP.Par.Train=0.\n\nPress any key to continue or cancel (CTRL+C).');        
    FP = Example1_FlatPlateAnalytical(FP,GP.Par.Lag); %FP model analytical
    GP.Par.Indx=sort(randperm(FP.Train.Samp,floor(FP.Train.Samp/GP.Par.SubSet))); %Subset of regressors for training  

    %Lift
    GP.Train(1).x=[FP.Train.alpha_h  FP.Train.alpha_a  FP.Train.alpha_h_d  FP.Train.alpha_a_d FP.Train.alpha_w  FP.Train.alpha_w_d]; %Organise training vector
    GP.Train(1).x_max=max(max(abs(GP.Train(1).x))); %Normalization value 
    GP.Train(1).x=GP.Train(1).x/GP.Train(1).x_max;  %Normalize input
 
    GP.Train(1).y_max=max(max(abs(FP.Train.CL+FP.Train.CLw))); %Normalize parameters 
    GP.Train(1).y=(FP.Train.CL+FP.Train.CLw)/GP.Train(1).y_max; %Normalization force
    GP.Train(1).noise=randn([FP.Train.Samp,1])/GP.Par.SNR*std(GP.Train(1).y); %Add noise for training.
    GP.Train(1).y=GP.Train(1).y+GP.Train(1).noise; %Add noise

    hyp=log(rand(8+3*GP.Par.Lag,1)); %Prior hyper parameters [Amp_kernel,Noise,alpha_h+Lag,alpha_a+Lag,alpha_h_d,alpha_a_d,alpha_w+Lag,alpha_w_d] in vectors.
    
    if GP.Par.Mean                   %Consider prior aerodynamic analytical mean
    GP.meanFlag(1).Mean = 1;         %Mean flag -> Consider aerodynamic prior
    GP.meanFlag(1).mfun=@Example1_Lift_LQS; %Mean aerodynamic model - Linear quasi steady model
    GP.meanFlag(1).x_max=GP.Train(1).x_max; %Maximum value required for normalisation of input for LQS
    GP.meanFlag(1).y_max=GP.Train(1).y_max; %Maximum value required for normalisation of output of LQS
    GP.meanFlag(1).Lag = GP.Par.Lag;        %Number of lag terms
    else
    GP.meanFlag(1).Mean = 0;         %Mean flag -> Do not consider aerodynamic prior
    end
    
    f=@(hyp)GP_Process_Opt(GP.Train(1).x(GP.Par.Indx,:),GP.Train(1).y(GP.Par.Indx),hyp,GP.Par.jitter,GP.Par.Noise,GP.meanFlag(1)); %Set up likelihood function to minimise
    GP.Train(1).hyp=fminunc(f,hyp,GP.Par.OptiControl1); % Minimize > get yperparameters

    %Moment
    GP.Train(2).x=GP.Train(1).x; %Same training input for moment as for lift
    GP.Train(2).x_max=GP.Train(1).x_max; %Normalization value > same as for lift

    GP.Train(2).y_max=max(max(abs(FP.Train.CM+FP.Train.CMw)));  %Normalize parameters  
    GP.Train(2).y=(FP.Train.CM+FP.Train.CMw)./GP.Train(2).y_max; %Normalization force 
    GP.Train(2).noise=randn([FP.Train.Samp,1])/GP.Par.SNR*std(GP.Train(2).y); %Add noise for training.
    GP.Train(2).y=GP.Train(2).y+GP.Train(2).noise; %Add noise

    if GP.Par.Mean                   %Consider prior aerodynamic analytical mean for moment (similar comments as for the moment)
    GP.meanFlag(2).Mean = 1;
    GP.meanFlag(2).mfun=@Example1_Moment_LQS;
    GP.meanFlag(2).x_max=GP.Train(2).x_max;
    GP.meanFlag(2).y_max=GP.Train(2).y_max;
    GP.meanFlag(2).Lag = GP.Par.Lag;
    else
    GP.meanFlag(2).Mean = 0; 
    end
    
    f=@(hyp)GP_Process_Opt(GP.Train(2).x(GP.Par.Indx,:),GP.Train(2).y(GP.Par.Indx),hyp,GP.Par.jitter,GP.Par.Noise,GP.meanFlag(2)); %Set up likelihood function to minimise
    GP.Train(2).hyp=fminunc(f,hyp,GP.Par.OptiControl1); % Minimize > get yperparameters
    save('Example1_FlatPlateAnalytical/Example1_Train','GP');%Save GP hyperparameters
else
    load('Example1_FlatPlateAnalytical/Example1_Train','GP');%Load GP hyperparameters 
end

%%  Prediction at Training Input
%Lift
GP.Pred(1).x=GP.Train(1).x; 
GP.Pred(1).y_targ=GP.Train(1).y-GP.Train(1).noise;
[GP.Pred(1).y,GP.Pred(1).y_var,~,GP.Train(1).L_Kern,~,~,~,~,~,GP.Pred(1).m]=GP_Process(GP.Train(1).x,GP.Train(1).y,GP.Pred(1).x,GP.Train(1).hyp,GP.Par.jitter,GP.Par.Noise,GP.meanFlag(1)); %Predictions

%Moment
GP.Pred(2).x=GP.Train(2).x;
GP.Pred(2).y_targ=GP.Train(2).y-GP.Train(2).noise;
[GP.Pred(2).y,GP.Pred(2).y_var,~,GP.Train(2).L_Kern,~,~,~,~,~,GP.Pred(2).m]=GP_Process(GP.Train(2).x,GP.Train(2).y,GP.Pred(2).x,GP.Train(2).hyp,GP.Par.jitter,GP.Par.Noise,GP.meanFlag(2)); %Predictions


%% Prediction at similar input as the training, only different random seed -> Consider moment only
rng(2); %Reproducibility
FP.Pred.Excitation='Train'; %Train 
FP = Example1_FlatPlateAnalytical(FP,GP.Par.Lag); %FP model analytical

GP.meanFlag(3) = GP.meanFlag(2);
GP.Train(3)=GP.Train(2);                          %Setup prediction model input & output
GP.Pred(3).x=[FP.Train.alpha_w*0  FP.Train.alpha_w*0  FP.Train.alpha_w_d*0  FP.Train.alpha_w_d*0 FP.Train.alpha_w  FP.Train.alpha_w_d]; %Initiate vector
GP.Pred(3).x=GP.Pred(3).x/GP.Train(3).x_max;      %Normalize input  
GP.Pred(3).y_targ=FP.Pred.CMw./GP.Train(3).y_max; %Normalize output
[GP.Pred(3).y,GP.Pred(3).y_var,GP.Pred(3).m]=GP_Process_Predict(GP.Train(3).x,GP.Train(3).y,GP.Pred(3).x,GP.Train(3).hyp,GP.Train(3).L_Kern,GP.meanFlag(3)); %Predictions
 
%% Prediction at Sinusoidal Motion Input
FP.Pred.Vr=[1:1:16];      %Prediction range for sinusoidal input motion - Vr range
FP.Pred.Alpha_Amp=1;      %Prediction angle amplitude 
FP.Pred.Cycl=ones(1,length(FP.Pred.Vr)).*6; %P Cycles - number of cycles for prediction
FP.Pred.Cycl(1)=12; % 12 cycles for Vr=1 to account for all lags.
StepsVr=cumsum([1 FP.Pred.Cycl.*FP.Pred.Vr./FP.Par.ds+1]); %Multi-step-ahead prediction

% Forced Vertical Displacements
FP.Pred.Excitation='SinH'; %SinH for sinusoidal input motion (heave)
FP = Example1_FlatPlateAnalytical(FP,GP.Par.Lag); %FP model analytical (for prediction)

%Lift
GP.meanFlag(4) = GP.meanFlag(1);
GP.Train(4)=GP.Train(1);
GP.Pred(4).x=[FP.Pred.alpha_h  FP.Pred.alpha_a  FP.Pred.alpha_h_d  FP.Pred.alpha_a_d FP.Pred.alpha_w FP.Pred.alpha_w_d];
GP.Pred(4).x=GP.Pred(4).x./GP.Train(4).x_max;              
GP.Pred(4).y_targ=FP.Pred.CL./GP.Train(4).y_max;
GP.Pred(4).y=GP.Pred(4).y_targ*0;
for i=1:length(FP.Pred.Vr)
 [GP.Pred(4).y(StepsVr(i):StepsVr(i+1)-1)]=GP_Process_Predict(GP.Train(4).x,GP.Train(4).y,GP.Pred(4).x((StepsVr(i):StepsVr(i+1)-1),:),GP.Train(4).hyp,GP.Train(4).L_Kern,GP.meanFlag(4)); %Predictions
end

%Moment
GP.meanFlag(5) = GP.meanFlag(2);
GP.Train(5)=GP.Train(2);
GP.Pred(5).x=GP.Pred(4).x;
GP.Pred(5).y_targ=FP.Pred.CM./GP.Train(5).y_max;
GP.Pred(5).y=GP.Pred(5).y_targ*0;
for i=1:length(FP.Pred.Vr)
 [GP.Pred(5).y(StepsVr(i):StepsVr(i+1)-1)]=GP_Process_Predict(GP.Train(5).x,GP.Train(5).y,GP.Pred(5).x((StepsVr(i):StepsVr(i+1)-1),:),GP.Train(5).hyp,GP.Train(5).L_Kern,GP.meanFlag(5)); %Predictions
end

% Forced Rotation
FP.Pred.Excitation='SinA'; %SinA for sinusoidal input motion (pitch)
FP = Example1_FlatPlateAnalytical(FP,GP.Par.Lag); %FP model analytical (for prediction)

%Lift
GP.meanFlag(6) = GP.meanFlag(1);
GP.Train(6)=GP.Train(1);
GP.Pred(6).x=[FP.Pred.alpha_h  FP.Pred.alpha_a  FP.Pred.alpha_h_d  FP.Pred.alpha_a_d FP.Pred.alpha_w FP.Pred.alpha_w_d];
GP.Pred(6).x=GP.Pred(6).x./GP.Train(6).x_max;              
GP.Pred(6).y_targ=FP.Pred.CL./GP.Train(6).y_max;
GP.Pred(6).y=GP.Pred(6).y_targ*0;
for i=1:length(FP.Pred.Vr)
 [GP.Pred(6).y(StepsVr(i):StepsVr(i+1)-1)]=GP_Process_Predict(GP.Train(6).x,GP.Train(6).y,GP.Pred(6).x((StepsVr(i):StepsVr(i+1)-1),:),GP.Train(6).hyp,GP.Train(6).L_Kern,GP.meanFlag(6)); %Predictions
end
 
%Moment
GP.meanFlag(7) = GP.meanFlag(2);
GP.Train(7)=GP.Train(2);
GP.Pred(7).x=GP.Pred(6).x;
GP.Pred(7).y_targ=FP.Pred.CM./GP.Train(7).y_max;
GP.Pred(7).y=GP.Pred(7).y_targ*0;
for i=1:length(FP.Pred.Vr)
 [GP.Pred(7).y(StepsVr(i):StepsVr(i+1)-1)]=GP_Process_Predict(GP.Train(7).x,GP.Train(7).y,GP.Pred(7).x((StepsVr(i):StepsVr(i+1)-1),:),GP.Train(7).hyp,GP.Train(7).L_Kern,GP.meanFlag(7)); %Predictions
end

%% Prediction at Static Input motion
FP.Pred.Excitation='Stat';             %Statatic wind coefficients
FP.Pred.StatAlpha=[-1 -.5 0 .5 1];     %Prediction angle amplitude
FP.Pred.StatTau=20;                    %Nondimensional time for prediction (it is later averaged)
Inc=(FP.Pred.StatTau./FP.Par.ds+1);    %Increment in the prediction vector

FP = Example1_FlatPlateAnalytical(FP,GP.Par.Lag); %FP model analytical (for prediction)

%Lift
GP.meanFlag(8) = GP.meanFlag(1);
GP.Train(8)=GP.Train(1);
GP.Pred(8).x=[FP.Pred.alpha_h  FP.Pred.alpha_a  FP.Pred.alpha_h_d  FP.Pred.alpha_a_d FP.Pred.alpha_h*0 FP.Pred.alpha_h_d*0];
GP.Pred(8).x=GP.Pred(8).x./GP.Train(8).x_max;              
GP.Pred(8).y_targ=FP.Pred.CL./GP.Train(8).y_max;
GP.Pred(8).y=GP.Pred(8).y_targ*0;
for i=1:length(FP.Pred.StatAlpha)
[GP.Pred(8).y(1+(i-1)*Inc:i*Inc)]=GP_Process_Predict(GP.Train(8).x,GP.Train(8).y,GP.Pred(8).x((1+(i-1)*Inc:i*Inc),:),GP.Train(8).hyp,GP.Train(8).L_Kern); %Predictions
end

%Moment
GP.meanFlag(9) = GP.meanFlag(2);
GP.Train(9)=GP.Train(2);
GP.Pred(9).x=[FP.Pred.alpha_h  FP.Pred.alpha_a  FP.Pred.alpha_h_d  FP.Pred.alpha_a_d FP.Pred.alpha_h*0 FP.Pred.alpha_h_d*0];
GP.Pred(9).x=GP.Pred(9).x./GP.Train(9).x_max;              
GP.Pred(9).y_targ=FP.Pred.CM./GP.Train(9).y_max;
GP.Pred(9).y=GP.Pred(9).y_targ*0;
for i=1:length(FP.Pred.StatAlpha)
[GP.Pred(9).y(1+(i-1)*Inc:i*Inc)]=GP_Process_Predict(GP.Train(9).x,GP.Train(9).y,GP.Pred(9).x((1+(i-1)*Inc:i*Inc),:),GP.Train(9).hyp,GP.Train(9).L_Kern); %Predictions
end

 %% Prediction at Sinusoidal Gust Input
FP.Pred.VrW=[1:2:30];       %Prediction range for sinusoidal input gust - Vr range
FP.Pred.Alpha_AmpW=0.5;     %Prediction angle amplitude 
FP.Pred.CyclW=ones(1,length(FP.Pred.VrW)).*6; %P Cycles - number of cycles for prediction
FP.Pred.CyclW(1)=12;         % 12 cycles for Vr=1 to account for all lags.
StepsVr=cumsum([1 FP.Pred.CyclW.*FP.Pred.VrW./FP.Par.ds+1]); %Multi-step-ahead prediction

FP.Pred.Excitation='SinW';  %SinW for sinusoidal input gust
FP = Example1_FlatPlateAnalytical(FP,GP.Par.Lag); %FP model analytical (for prediction)
%Lift
GP.meanFlag(10) = GP.meanFlag(1);
GP.Train(10)=GP.Train(1);
GP.Pred(10).x=[FP.Pred.alpha_w*0  FP.Pred.alpha_w*0  FP.Pred.alpha_w_d*0  FP.Pred.alpha_w_d*0 FP.Pred.alpha_w  FP.Pred.alpha_w_d];
GP.Pred(10).x=GP.Pred(10).x./GP.Train(10).x_max;              
GP.Pred(10).y_targ=FP.Pred.CLw./GP.Train(10).y_max;
GP.Pred(10).y=GP.Pred(10).y_targ*0;
for i=1:length(FP.Pred.VrW)
  [GP.Pred(10).y(StepsVr(i):StepsVr(i+1)-1),~,m(StepsVr(i):StepsVr(i+1)-1)]=GP_Process_Predict(GP.Train(10).x,GP.Train(10).y,GP.Pred(10).x((StepsVr(i):StepsVr(i+1)-1),:),GP.Train(10).hyp,GP.Train(10).L_Kern,GP.meanFlag(10)); %Predictions
end
 
%Moment
GP.meanFlag(11) = GP.meanFlag(2);
GP.Train(11)=GP.Train(2);
GP.Pred(11).x=[FP.Pred.alpha_w*0  FP.Pred.alpha_w*0  FP.Pred.alpha_w_d*0  FP.Pred.alpha_w_d*0 FP.Pred.alpha_w  FP.Pred.alpha_w_d];
GP.Pred(11).x=GP.Pred(11).x./GP.Train(11).x_max;              
GP.Pred(11).y_targ=FP.Pred.CMw./GP.Train(11).y_max;
GP.Pred(11).y=GP.Pred(11).y_targ*0;
for i=1:length(FP.Pred.VrW)
  [GP.Pred(11).y(StepsVr(i):StepsVr(i+1)-1),~,m(StepsVr(i):StepsVr(i+1)-1)]=GP_Process_Predict(GP.Train(11).x,GP.Train(11).y,GP.Pred(11).x((StepsVr(i):StepsVr(i+1)-1),:),GP.Train(11).hyp,GP.Train(11).L_Kern,GP.meanFlag(11)); %Predictions
end

% save('Example1_FlatPlateAnalytical/Ex1c_Out','GP','FP','-v7.3');
%% Plots
% load('Example1_FlatPlateAnalytical/Ex1c_Out','GP','FP');
Example1c_Plots(GP,FP); 
