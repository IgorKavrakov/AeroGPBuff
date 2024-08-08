function [Mean_Pred,Var,m_pred,lik] = GP_Process_Predict(x,y,x_predict,hyp,L_Kern,meanFlag)
% Prediction for a GP - similar as GP_Process.m, excluding the computation of kernel and likelihood.

% By Igor Kavrakov

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

%%
if nargin<6||meanFlag.Mean==-1
   m=y.*0;m_pred=0;
else
   m=meanFlag.mfun(x,hyp,meanFlag); 
   m_pred=meanFlag.mfun(x_predict,hyp,meanFlag);  
end

Kern_predict_r=Kern_Exp(x_predict,x,hyp);

alpha=L_Kern'\(L_Kern\(y-m));
Mean_Pred=Kern_predict_r*alpha+m_pred;

if nargout>1
Kern_predict=Kern_Exp(x_predict,x_predict,hyp); %Get the kernel at prediction points    
v=L_Kern\Kern_predict_r';       %Predictive varianse
Var=Kern_predict-v'*v;          %Predictive variances
elseif nargout>3   
lik=-1/2*(y-m)'*alpha-sum(log(diag(L_Kern)))-length(y)/2*log(2*pi); %Likelihood   
end

end

