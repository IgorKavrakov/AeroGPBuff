function [Mean_Pred,m_pred] = GP_Process_Predict_Mean(x,~,x_predict,hyp,alpha,meanFlag)
% Prediction for a the mean function for a GP

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
if nargin<6||meanFlag.Mean==-1
   m_pred=0;
else
   m_pred=meanFlag.mfun(x_predict,hyp,meanFlag);  
end

Kern_predict_r=Kern_Exp(x_predict,x,hyp);

Mean_Pred=Kern_predict_r*alpha+m_pred;

end

