% Find the optimal Sigma and Volumn factor
% Date: Oct 8th, 2019
% Author: Yang LI
close all
clear

%% Initialization
excavator_data

% parameters = zeros(3,7);

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
depH = depH - 0.5*min(depH,[],"all");

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

%% Error analysis
error_before = immse(depH, ModelH);

%% Optimization
theta0 = [2000,2000,3];
lb = [1000,1400,2.5];
ub = [2500,3000,3.5];

fun = @(theta)immse(depH, function_input_2d(X,Y,c,theta(3)*V,[theta(1),0;0,theta(2)],the,xf,yr,yl));

options = optimset('Display','iter','PlotFcns',@optimplotfval);

tic
% theta = fmincon(fun,theta0,[],[],[],[],lb,ub);
theta = fminsearch(fun,theta0);
toc

parameters(:,i) = theta';

%% result
lambdaX = theta(1); lambdaY = theta(2); kV = theta(3);

ModelH = function_input_2d(X,Y,c,kV*V,Sigma,the,xf,yr,yl);

error_before
error_after = immse(depH, ModelH)

% figure
% mesh(X,Y,deltaH)
% xlabel('x[mm]')
% ylabel('y[mm]')
% zlabel('h[mm]')
% xlim([0 170])
% zlim([0 40])

%% error model
% H_error = depH - ModelH;
% 
% figure
% mesh(X,Y,H_error)
% xlabel('x[mm]')
% ylabel('y[mm]')
% zlabel('h[mm]')
% xlim([0 170])
% zlim([-10 10])

% input('Next data?[Enter]')
end

%% summary
maybeGoodResult = mean(parameters,2)
% % resultant parameters are obtained from 7 data
