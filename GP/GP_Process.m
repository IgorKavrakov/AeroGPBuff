function [Mean_Pred,Var,lik,L_Kern,Kern,Kern_predict,Kern_predict_r,alpha,v,m_pred] = GP_Process(x,y,x_predict,hyp,jitter,noise,meanFlag)
%This is implementation of GP regression based on Rasmussen&Willams Gaussian Processes for Machine Learning (2006) (Algorithm 2.1)

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

%%
if noise~=-1
    noise=log(0);
else
    noise=hyp(end);    
end

if nargin<7||meanFlag.Mean==-1
   m=y.*0;m_pred=0;
else
   m=meanFlag.mfun(x,hyp,meanFlag); 
   m_pred=meanFlag.mfun(x_predict,hyp,meanFlag); 
end

Kern=Kern_Exp(x,x,hyp); %Get the kernel at x points
Kern_predict=Kern_Exp(x_predict,x_predict,hyp); %Get the kernel at prediction points
Kern_predict_r=Kern_Exp(x_predict,x,hyp); %Get kernel at cross points
Noise=eye(size(x,1)).*exp(noise)+eye(size(x,1)).*jitter;

L_Kern=chol(Kern+Noise,'lower'); %Compute cholesky (for inversion - faster). Check Rasmussen
alpha=L_Kern'\(L_Kern\(y-m));       %Predictive mean
Mean_Pred=Kern_predict_r*alpha+m_pred; %Predictive the mean
v=L_Kern\Kern_predict_r';       %Predictive varianse
Var=Kern_predict-v'*v;          %Predictive variances
lik=-1/2*(y-m)'*alpha-sum(log(diag(L_Kern)))-length(y)/2*log(2*pi); %Likelihood


end

