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
% u = [xc1,yc1,vol1, xc2,yc2,vol2, xc3,yc3,vol3];
% u0 = [85,0,1, 85,0,1, 85,0,1];
% u0 = [42.5,0,1, 85,0,1, 3*42.5,0,1];
% u0 = [42.5,0,1, 85,0,1, 3*42.5,0,1];
u0 = [42.5000000000000,-6.86869124880580,1,87.5499824485590,-0.674760706942117,1,121.549997376126,2.04452083448939,0.999999966438559];

% boundaries:
umin = [42.5,-40,0, 42.5,-40,0, 42.5,-40,0];
umax = [3*42.5,40,1, 3*42.5,40,1, 3*42.5,40,1];

% constraints;
A = []; b = []; Aeq = []; beq = [];

%% 3. Optimization
tic
% options = optimoptions(@fmincon,'MaxFunctionEvaluations',6000,'Display','iter','Algorithm','sqp');
% u_opt = fmincon(@(u)objfun_3times_loading(u,H0,R,Nominal_model,X_data,Y_data,X,Y,hyp_sparseGP,U),...
%     u0,A,b,Aeq,beq,umin,umax,[],options);

options = optimset('Display','iter','PlotFcns',@optimplotfval);
u_opt = fminsearch(@(u)objfun_3times_search(u,H0,R,Nominal_model,X_data,Y_data,X,Y,hyp_sparseGP,U),...
    u0,options);
toc

%% 4. Show the result
% Loading soil to the end
H_last = H0;
for i = 1:3
    H_after = gp_predict(H_last,u_opt(3*i-2),u_opt(3*i-1),u_opt(3*i),U,X,Y,...
        X_data,Y_data,hyp_sparseGP, Nominal_model);
    H_last = H_after;
end
figure
mesh(H_after)
figure
subplot(1,2,1); mesh(H_after);zlim([-50 30]);
subplot(1,2,2); mesh(R);zlim([-50 30]);

% 5. save
% save("optimization results\"+date+"-u3_opt",'u_opt')