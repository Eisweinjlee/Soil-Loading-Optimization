%% Episodic loading planning
% Loading planning optimization with GP-based soil model
% Author: Li, Yang
% Date: Mar 9th, 2020

% close all
% clear

nominal_model_flag = 1; 
% 0: 2D Gaussian pdf, 1: Unit-symmetric Gaussian pdf

%% 1. Parameters, dataset and model
run("data\excavator_data.m") % load parameters

data_date = "09-Mar-2020-";
load("dataset\"+data_date+"local_dataset.mat") % load the dataset
% Included: X_data, Y_data, X, Y
load("trained GPs\"+data_date+"model.mat") % load trained parameters
% Included: hyp_sparseGP

% two Nominal models
if nominal_model_flag == 0 % 0: 2D Gaussian pdf
    lambdaX = 2.302227968708186e+03; lambdaY = 1.631169900899356e+03;
    kV = 3.445999018892074; Xmove = 2.1830e-07; Ymove = 0.207118001857057;
    Nominal_model = [nominal_model_flag,Xmove,Ymove,kV,lambdaX,lambdaY,the,xf,yr,yl];
% function_input_2d(X,Y,[Xc(i);Yc(i)]-[Xmove;Ymove],kV*Vol_data(i)*Vol,Sigma,the,xf,yr,yl);

elseif nominal_model_flag == 1 % 1: Unit-symmetric Gaussian pdf
    lambdaX = 1.209827074460249e+03; lambdaY = 1.857790332589524e+03;
    kV = 3.9410; Xmove = 0.0253; Ymove = 0.1487; alpha = 0.6147;
    Nominal_model = [nominal_model_flag,Xmove,Ymove,kV,lambdaX,lambdaY,alpha];
% unit_sym_input_2d(X,Y,[Xc(i);Yc(i)]-[Xmove,Ymove],kV*Vol_data(i)*Vol,Sigma,alpha);
end

load("data\pca data\U_matrix_PCA.mat") % load U matrix for PCA compute

%% 2. Optimization parameters
% inputs: each loaing has (Xc,Yc,Vol)
% (format) u = [xc1,yc1,vol1, xc2,yc2,vol2, xc3,yc3,vol3, .....];
% u0 = [42.5,0,1, 85,0,1, 3*42.5,0,1];
% u0 = [42.51,-7.483,1, 87.262,-0.851,1, 121.55,7.621,0.973]; % search result

% Initial inputs
% u0 = [42.51,0,1, 87.262,0,0.5, 100,0,0.6, 125,0,0.8];
% u0 = [53.5499989365922,-6.78790062206207,0.863609289821117,67.1499930688198,6.51546386583886,0.719505660957496,86.5122258799381,-3.55100792616697,0.797932860498399,125.404269206472,9.19056480218624,0.832308210767386];
u0 = [53.5493935897970,-3.87599608955972,0.831205387328661,67.1450471157161,8.72476005996204,0.735134898326868,85.8510758456105,-3.62645944766534,0.770890789760168,126.670983679879,9.36170179718247,0.813616886168327];

% constraints;
umin = []; umax = [];
loading_times = length(u0)/3;
for i = 1:loading_times
    umin = [umin,42.5,-40,0]; % [Xc_min, Yc_min, Vol_min]
    umax = [umax,127.5,40,1]; % [Xc_max, Yc_max, Vol_max]
end
A = []; b = []; Aeq = []; beq = [];

%% 3. Optimization
tic
options = optimoptions(@fmincon,'MaxFunctionEvaluations',6000,'Display','iter','Algorithm','sqp');
u_opt = fmincon(@(u)objfun_loading(u,H0,R,Nominal_model,X_data,Y_data,X,Y,hyp_sparseGP,U),...
    u0,A,b,Aeq,beq,umin,umax,[],options);

% options = optimset('Display','iter','PlotFcns',@optimplotfval,'MaxFunEvals',6000);
% u_opt = fminsearch(@(u)objfun_search(u,H0,R,Nominal_model,X_data,Y_data,X,Y,hyp_sparseGP,U),...
%     u0,options);

y_iszero = 0;

% % (Optional) Y is set to be zero to guarantee the symmetricity 
% u0 = [43,1, 95.2,0.7, 123.25,0.49, 124.15,0.72];
% umin = [42.5,0, 42.5,0, 42.5,0, 42.5,0];
% umax = [127.5,1, 127.5,1, 127.5,1, 127.5,1];
% options = optimoptions(@fmincon,'MaxFunctionEvaluations',6000,'Display','iter','Algorithm','sqp');
% u_opt = fmincon(@(u)objfun_loading_symAtY(u,H0,R,Nominal_model,X_data,Y_data,X,Y,hyp_sparseGP,U),...
%     u0,A,b,Aeq,beq,umin,umax,[],options);y_ iszero = 1;

toc
%% 4. Show the result
% Loading soil to the end
H_last = H0;
if y_iszero == 1
    for i = 1:loading_times
    H_after = gp_predict(H_last,u_opt(2*i-1),0,u_opt(2*i),U,X,Y,...
        X_data,Y_data,hyp_sparseGP, Nominal_model);
    H_last = H_after;
    end
    
elseif y_iszero == 0
    for i = 1:loading_times
        H_after = gp_predict(H_last,u_opt(3*i-2),u_opt(3*i-1),u_opt(3*i),U,X,Y,...
            X_data,Y_data,hyp_sparseGP, Nominal_model);
        H_last = H_after;
    end   
end

figure
mesh(X,Y,H_after)
figure
subplot(1,2,1); mesh(X,Y,H_after);xlim([0 170]);ylim([-80 80]);zlim([-50 30]);
subplot(1,2,2); mesh(X,Y,R);xlim([0 170]);ylim([-80 80]);zlim([-50 30]);


% 5. save
% save("optimization results\"+date+"-u4_opt",'u_opt')