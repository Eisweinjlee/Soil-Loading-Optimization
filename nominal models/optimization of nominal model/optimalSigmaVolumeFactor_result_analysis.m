% Try the optimal Sigma and Volumn factor
% Date: Oct 8th, 2019
% Author: Yang LI
close all
clear

%% Initialization
excavator_data

parameters = zeros(3,7);

for i = 1:7
filename = "data_centerdep/data"+ i +".mat";
load(filename) % data1.mat ~ data7.mat
H = dep; clear dep

% figure
% mesh(X,Y,H)
% xlabel('x[mm]')
% ylabel('y[mm]')
% zlabel('h[mm]')
% zlim([-50 40])

depH = H - H0;
depH = depH - 0.5 * min(depH,[],"all");

% figure
% mesh(X,Y,depH)
% xlabel('x[mm]')
% ylabel('y[mm]')
% zlabel('h[mm]')
% xlim([0 170])
% zlim([0 40])

[row,col] = find(depH == max(depH,[],'all')); % the peak

%% optimization program initialization
lambdaX = 1709.10978; lambdaY = 2274.09987; kV = 2.9126;
% our best parameter for Gaussian distribution

Sigma = [lambdaX,0; 0,lambdaY];
depx = X(1,col); depy = Y(row,1); c = [depx, depy]; 
the = atan((depy-Pe(2))/(depx-Pe(1)));
V = sum(depH ,'all'); % 1.8e+5 is the actual volume

ModelH = function_input_2d(X,Y,c,kV*V,Sigma,the,xf,yr,yl);

% figure
% mesh(X,Y,deltaH)
% xlabel('x[mm]')
% ylabel('y[mm]')
% zlabel('h[mm]')
% xlim([0 170])
% zlim([0 40])

%% Error model analysis
error = immse(depH, ModelH);

H_error(:,:,i) = depH - ModelH;

titlename = "Data No."+ i +", MSE of error is "+ error +".";

figure
mesh(X,Y,H_error(:,:,i))
xlabel('x[mm]')
ylabel('y[mm]')
zlabel('h[mm]')
xlim([0 170])
zlim([-10 10])
title(titlename)
end