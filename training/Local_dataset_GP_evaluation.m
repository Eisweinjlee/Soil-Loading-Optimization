%% Local_dataset_GP_evaluation
% The trained GP model for error distribution is evaluated.
% Author: Li, Yang
% Date: Jan 15th, 2020
% Modification: Mar 7th, 2020
close all
clear

nominal_model_flag = 1; % 0: 2D Gaussian pdf, 1: Unit-symmetric Gaussian pdf

%% 1. Load the model and specify parameters

data_date = "09-Mar-2020-";

% load the dataset: X_data, Y_data, X, Y
load("dataset\"+data_date+"local_dataset.mat")
% load trained parameters: hyp_sparseGP
load("trained GPs\"+data_date+"model.mat") % 800 virtual data points
    % load("trained GPs\"+data_date+"model.mat")

% Specify the mean, cov, likelihood
meanfunc = [];                      % empty: don't use a mean function
covfunc = {@covSEard};              % ARD SE
likfunc = {@likGauss};              % Gaussian likelihood
infmethod = @infGaussLik;           % inference with Guassian Likelihood

% inducing points
xu = hyp_sparseGP.xu; cov = {'apxSparse', covfunc, xu};
inff = @(varargin) infmethod(varargin{:},struct('s', 0));
% VFE, opt.s = 0; SPEP, 0 <opt.s < 1; FITC, opt.s = 1

%% 2. Load the data for evaluation
m = 94; n = 100;

% loading centers
load("dataset\"+data_date+"loading_center.mat")
% Volumes
load("dataset\"+data_date+"loading_volume.mat")
% H_data: full_dataset
load("dataset\"+data_date+"full_dataset.mat")
number = length(H_data(1,1,:));

% shape before loading: PCs & mean value
load("dataset\"+data_date+"pca_mean.mat")
H_last = zeros(94,100,number); 
for i = 1:97
    H_last(:,:,i) = H0; 
end
number = 97;
for i = 1:49
    H_last(:,:,number+1) = H0;
    H_last(:,:,number+2) = H_data(:,:,number+1);
    number = number + 2;
end
for i = 1:10
    H_last(:,:,number+1) = H0;
    H_last(:,:,number+2) = H_data(:,:,number+1);
    H_last(:,:,number+3) = H_data(:,:,number+2);
    number = number + 3;
end


%% 3. The prediction
err = zeros(number,1);
err_nosigmoid = zeros(number,1);
err_noGP = zeros(number,1);

if nominal_model_flag == 0
    lambdaX = 2.302227968708186e+03; lambdaY = 1.631169900899356e+03;
    kV = 3.445999018892074; Xmove = 2.1830e-07; Ymove = 0.207118001857057;
    Sigma = [lambdaX,0;0,lambdaY];
elseif nominal_model_flag == 1
    lambdaX = 1.209827074460249e+03; lambdaY = 1.857790332589524e+03;
    kV = 3.9410; Xmove = 0.0253; Ymove = 0.1487; alpha = 0.6147;
    Sigma = [lambdaX,0;0,lambdaY];
end

for i = 1:number % start evaluation
if mod(i,5) == 0
    disp("Evaluation process "+i+"/"+number)
end

% for i = 185

% normalized center: Xc=[0,1], Yc=[-1,1]
Xc_nor = Xc(i)/170; Yc_nor = Yc(i)/80;

% pca_mean normalization
pca_mean_vec = pca_mean(i,:);
pca_mean_vec = ones(9400,1)*[pca_mean_vec(1:4)/800, pca_mean_vec(4+1)/40];

% test data
X_test = [Xc_nor*ones(9400,1), Yc_nor*zeros(9400,1), ...
    (X(:)-Xc(i))/170, (Y(:)-Yc(i))/80, Vol_data(i)*ones(9400,1), pca_mean_vec];

% prediction
tic
[ymu,ys2] = gp(hyp_sparseGP, inff, meanfunc, cov, likfunc,...
    X_data, Y_data, X_test);
toc

% reshape the result
m = 94; n = 100;
H_error_pred = reshape(ymu,[m,n]);
CovDist = reshape(ys2,[m,n]);

% % Plot: the mean prediction
% figure; mesh(X, Y, H_error_pred)
% title("$X_c=$"+Xc*170+", $Y_c=$"+Yc*80,'Interpreter','latex')
% 
% % Plot: covariance analysis
% s_upper = ymu + 2*sqrt(ys2); s_lower = ymu - 2*sqrt(ys2);
% s_upper = reshape(s_upper,[m,n]); s_lower= reshape(s_lower,[m,n]);
% figure; hold on;
% sur1 = surf(X,Y,H_error); sur1.EdgeColor = 'flat';
% sur2 = surf(X,Y,s_upper,'FaceAlpha',0.3); sur2.EdgeColor = 'none';
% sur3 = surf(X,Y,s_lower,'FaceAlpha',0.3); sur3.EdgeColor = 'none';
% view([-23,20]); hold off;
% 
% figure; mesh(X,Y,CovDist);

%% 3. Relative covariance modification

% find the data close to center X(1,:), Y(:,1)
disX = (X(1,:)-Xc(i)).^2; disY = (Y(:,1)-Yc(i))'.^2;
[minX_v,minX_p] = min(disX); [minY_v,minY_p] =  min(disY);
data_cen = [X(1,minX_p);Y(minY_p,1)];

% % Plot: the loading center and data center
% figure;hold on;
% plot(X,Y,'LineStyle','none','Marker','.','Color',[218,165,32]/255)
% plot(Xc*170,Yc*80,'LineStyle','none','Marker','+','Color',[199,21,133]/255,'LineWidth',1.5)
% plot(data_cen(1),data_cen(2),'LineStyle','none','Marker','+','Color',[0,0,128]/255,'LineWidth',1.5)
% hold off;

% calculate all relative covariance
cov_cen = CovDist(minY_p,minX_p); % the covariance of center
w = cov_cen ./ CovDist;

% Sigmoid function
a = 1; b = 18.1; % https://www.desmos.com/calculator/kn9tpwdan5
% a = -0.6; b = 5.4;
gw = 1./(1 + exp(a + -b.*w)); % sigmoid
% q = 0:0.01:max(w,[],'all'); % plot to see the sigmoid
% gw = 1./(1 + exp(a + -b.*q)); % sigmoid
% figure; plot(q,gw); grid on

% modify the mean predition through sigmoid
H_modified = gw .* H_error_pred;
H_no_modified = H_error_pred;

% % Plot: before modification and after
% figure; subplot(1,2,1); mesh(X,Y,H_error_pred); zlim([-20 20]); title("before");
% subplot(1,2,2); mesh(X,Y,H_modified); zlim([-20 20]);title("after");

%% 4. Add to the nominal Gaussian pdf model
excavator_data;
Vol = 7.5889e+04;

if nominal_model_flag == 0
    H_nominal = function_input_2d(X,Y,[Xc(i);Yc(i)]-[Xmove;Ymove],kV*Vol_data(i)*Vol,Sigma,the,xf,yr,yl);
elseif nominal_model_flag == 1
    H_nominal = unit_sym_input_2d(X,Y,[Xc(i);Yc(i)]-[Xmove,Ymove],kV*Vol_data(i)*Vol,Sigma,alpha);
end

H_noGP = H_last(:,:,i) + H_nominal;

H_combination = H_noGP + H_modified;
H_nosigmoid = H_noGP + H_no_modified;

%% 5. Evaluate the performance
err(i) = sqrt(immse(H_combination, H_data(:,:,i)));
err_nosigmoid(i) = sqrt(immse(H_nosigmoid, H_data(:,:,i)));
err_noGP(i) = sqrt(immse(H_noGP, H_data(:,:,i)));

% figure; 
% subplot(1,3,1); mesh(X,Y,H_noGP); 
% xlim([0 170]); zlim([-50 40]); title("Nominal");zlabel("h[mm]")
% subplot(1,3,2); mesh(X,Y,H_combination); 
% xlim([0 170]); zlim([-50 40]); title("Nominal+GP");
% subplot(1,3,3); mesh(X,Y,H_data(:,:,i)); zlim([-20 20]);
% xlim([0 170]); zlim([-50 40]); title("Real, nMSE= " + err(i));
% % close all
% 
% figure; mesh(X,Y,H_combination); 
% xlim([0 170]); zlim([-50 15]); view([48 36]);
% figure; mesh(X,Y,H_data(:,:,i));
% xlim([0 170]); zlim([-50 15]); view([48 36]);

end

%% 6. plot the normalize MSE

mean(err)
figure; hold on; grid on;
plot(err,'LineStyle','none','Marker','x','Color','[0, 0, 0]','LineWidth',1.2);
plot([0 number],[mean(err) mean(err)],'--','Linewidth',1.2)
xlim([0 number]);  % ylim([0 15]);
ylabel("rMSE[mm]");

figure; hold on; grid on;
plot(err_noGP,'LineStyle','none','Marker','x','MarkerSize',11,'Color',[220,20,60]/255,'LineWidth',1.5);
plot(err_nosigmoid,'LineStyle','none','Marker','x','MarkerSize',11,'Color',[218,165,32]/255,'LineWidth',1.5);
plot(err,'LineStyle','none','Marker','x','MarkerSize',11,'Color',[0,0,139]/255,'LineWidth',1.5);
legend("Gaussian pdf","no mapping","result")
xlim([0 number]);  ylim([0 10]);
ylabel("rMSE[mm]");