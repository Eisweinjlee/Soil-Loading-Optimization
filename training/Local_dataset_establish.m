%% Establishing local dataset for loading data
% with the "data1_positions", the "data2_volumes"
% and the "data3_beforeLoading"
% Author: Li, Yang
% Date: Dec 3rd, 2019
% Modification: Mar 7th, 2020
close all
clear

run('data\excavator_data.m')
docName1 = "data\data1_positions\";
docName2 = "data\data2_volumes\";
docName3 = "data\data3_beforeLoading\";
save_doc = "dataset\";

nominal_model_flag = 1; % 0: 2D Gaussian pdf, 1: Unit-symmetric Gaussian pdf
rank = 4; % PCA for h_{i-1} (enough to use r = 4 in our project)

%% 1. Loading the data (changed for different data)
load(docName1 + "00.mat",'dep') % load the empty vessel data
H0 = dep; [m,n] = size(dep);
H_data = zeros(m,n,61);

%% 1.1 load docName1 - Position difference
number = 0;
for i = 1:9
    for j = 1:5
        number = number + 1;
        filename = docName1 + num2str(10*i+j)+".mat";
        load(filename)
        H_data(:,:,number) = dep;
    end
end
for i = 1:16
    number = number + 1;
    filename = docName1 + num2str(100+i)+".mat";
    load(filename)
    H_data(:,:,number) = dep;
end
% 61 profiles

%% 1.2 load docName2 - Volume difference
for i = 1:9
    for j = 1:4
        number = number + 1;
        filename = docName2 + num2str(10*i+j)+".mat";
        load(filename)
        H_data(:,:,number) = dep;
    end
end
% 36 profiles, 97 profiles included docName1

%% 1.3 load docName3 - H{i-1} dataset
load(docName3 + "00"); H02 = dep;

% b10 ~ b61
for i = 1:6
    filename = docName3 + "b" + num2str(10*i)+".mat"; load(filename);
    H_data(:,:,number + 1) = dep;
    filename = docName3 + "b" + num2str(10*i+1)+".mat"; load(filename);
    H_data(:,:,number + 2) = dep;
%     figure;subplot(1,2,1);mesh(H_data(:,:,number + 1))
%     subplot(1,2,2);mesh(H_data(:,:,number + 2))
    number = number + 2;
end
% bb10 ~ bb81
for i = 1:8
    filename = docName3 + "bb" + num2str(10*i)+".mat"; load(filename);
    H_data(:,:,number + 1) = dep;
    filename = docName3 + "bb" + num2str(10*i+1)+".mat"; load(filename);
    H_data(:,:,number + 2) = dep;
%     figure;subplot(1,2,1);mesh(H_data(:,:,number + 1))
%     subplot(1,2,2);mesh(H_data(:,:,number + 2))
    number = number + 2;
end
% c10 ~ c91
for i = 1:9
    filename = docName3 + "c" + num2str(10*i)+".mat"; load(filename);
    H_data(:,:,number + 1) = dep;
    filename = docName3 + "c" + num2str(10*i+1)+".mat"; load(filename);
    H_data(:,:,number + 2) = dep;
%     figure;subplot(1,2,1);mesh(H_data(:,:,number + 1))
%     subplot(1,2,2);mesh(H_data(:,:,number + 2))
    number = number + 2;
end
% l10 ~ l61
for i = 1:6
    filename = docName3 + "l" + num2str(10*i)+".mat"; load(filename);
    H_data(:,:,number + 1) = dep;
    filename = docName3 + "l" + num2str(10*i+1)+".mat"; load(filename);
    H_data(:,:,number + 2) = dep;
%     figure;subplot(1,2,1);mesh(H_data(:,:,number + 1))
%     subplot(1,2,2);mesh(H_data(:,:,number + 2))
    number = number + 2;
end
% r10 ~ r61
for i = 1:6
    filename = docName3 + "r" + num2str(10*i)+".mat"; load(filename);
    H_data(:,:,number + 1) = dep;
    filename = docName3 + "r" + num2str(10*i+1)+".mat"; load(filename);
    H_data(:,:,number + 2) = dep;
%     figure;subplot(1,2,1);mesh(H_data(:,:,number + 1))
%     subplot(1,2,2);mesh(H_data(:,:,number + 2))
    number = number + 2;
end
% t10 ~ t61
for i = 1:6
    filename = docName3 + "t" + num2str(10*i)+".mat"; load(filename);
    H_data(:,:,number + 1) = dep;
    filename = docName3 + "t" + num2str(10*i+1)+".mat"; load(filename);
    H_data(:,:,number + 2) = dep;
%     figure;subplot(1,2,1);mesh(H_data(:,:,number + 1))
%     subplot(1,2,2);mesh(H_data(:,:,number + 2))
    number = number + 2;
end
% tt10 ~ tt81
for i = 1:8
    filename = docName3 + "tt" + num2str(10*i)+".mat"; load(filename);
    H_data(:,:,number + 1) = dep;
    filename = docName3 + "tt" + num2str(10*i+1)+".mat"; load(filename);
    H_data(:,:,number + 2) = dep;
%     figure;subplot(1,2,1);mesh(H_data(:,:,number + 1))
%     subplot(1,2,2);mesh(H_data(:,:,number + 2))
    number = number + 2;
end
% 11 ~ 103
for i = 1:10
    filename = docName3 + num2str(10*i+1)+".mat"; load(filename);
    H_data(:,:,number + 1) = dep;
    filename = docName3 + num2str(10*i+2)+".mat"; load(filename);
    H_data(:,:,number + 2) = dep;
    filename = docName3 + num2str(10*i+3)+".mat"; load(filename);
    H_data(:,:,number + 3) = dep;
%     figure;subplot(1,3,1);mesh(H_data(:,:,number + 1))
%     subplot(1,3,2);mesh(H_data(:,:,number + 2))
%     subplot(1,3,3);mesh(H_data(:,:,number + 3))
    number = number + 3;
end
% 128 profiles, 225 profiles included docName1

filename = date + "-full_dataset.mat";
save(save_doc + filename,'H_data')

%% 2.1 The loading center data
% docName1
Xc = [25*ones(5,1);50*ones(5,1);75*ones(15,1);50*ones(10,1);25*ones(10,1)];
Yc = [zeros(15,1);40*ones(5,1);-40*ones(10,1);40*ones(10,1);-40*ones(5,1)];
Xc = [Xc;62.5;37.5;25;25;37.5;62.5;75;75;50;50;62.5;37.5;37.5;37.5;62.5;62.5];
Yc = [Yc;40;40;20;-20;-40;-40;-20;20;20;-20;0;0;-20;20;20;-20];

% docName2
Xc = [Xc;25*ones(4,1);50*ones(4,1);75*ones(12,1);50*ones(8,1);25*ones(8,1)];
Yc = [Yc;zeros(12,1);40*ones(4,1);-40*ones(8,1);40*ones(8,1);-40*ones(4,1)];

% docName3: b, bb, c, l, r, t, tt, 1x~10x
Xc = [Xc; 50;75; 75;75; 75;75; 75;75; 50;75; 50;75];
Yc = [Yc; 0;0; 0;0; 40;0; -40;0; -40;0; 40;0 ]; %b
Xc = [Xc; 50;75; 75;75; 75;75; 50;75; 50;75; 75;75; 75;75; 50;75];
Yc = [Yc; 0;40; 0;40; 40;40; 40;40; 0;-40; 0;-40; -40;-40; -40;-40];%bb
Xc = [Xc; 25;50; 50;50; 75;50; 75;50; 75;50; 50;50; 50;50; 25;50; 25;50];
Yc = [Yc; 0;0; 0;0; 0;0; 40;0; -40;0; -40;0; 40;0; 40;0; -40;0];%c
Xc = [Xc; 25;50; 50;50; 75;50; 75;50; 50;50; 25;50];
Yc = [Yc; 0;-40; 0;-40; 0;-40; -40;-40; -40;-40; -40;-40];%l
Xc = [Xc; 25;50; 50;50; 75;50; 75;50; 50;50; 25;50];
Yc = [Yc; 0;40; 0;40; 0;40; 40;40; 40;40; 40;40];%r
Xc = [Xc; 25;25; 50;25; 50;25; 50;25; 25;25; 25;25];
Yc = [Yc; 0;0; 0;0; 40;0; -40;0; -40;0; 40;0];%t
Xc = [Xc; 25;25; 50;25; 50;25; 25;25; 25;25; 50;25; 50;25; 25;25];
Yc = [Yc; 0;40; 0;40; 40;40; 40;40; 0;-40; 0;-40; -40;-40; -40;-40];%tt
Xc = [Xc; 50;50;25; 25;25;50; 25;50;75; 75;50;25; 25;75;50; 75;25;50;
    25;75;50; 25;75;50; 75;25;50; 75;25;50];
Yc = [Yc; 0;0;0; 0;0;0; 0;0;0; 0;0;0; 0;0;0; 0;0;0; 40;-40;0; -40;40;0;
    40;40;0; -40;-40;0];%1x~10x
Xc = Xc*1.7; % Xc\in[0,170]

filename = date + "-loading_center.mat";
save(save_doc + filename,'Xc','Yc')

%% 2.2 The volume data
Vol_data = ones(61,1); % docName1
for i = 1:9 
    Vol_data = [Vol_data;0.25;0.50;0.50;0.75]; % docName2
end
Vol_data = [Vol_data; ones(128,1)]; % docName3

filename = date + "-loading_volume.mat";
save(save_doc + filename,'Vol_data')

%% 2.3 The shape before loading - PCA + mean
load('Soil Loading Optimization\data\pca data\U_matrix_PCA.mat')
pca_mean = zeros(number,rank+1);

for i = 1:97
    pca_mean(i,:) = PCA_pc([Xc(i);Yc(i)],H0,rank,U);
end % docName1 & docName2

% docName3:b,bb,c,l,r,t,tt,1x~10x
number = 97;
for i = 1:6
    pca_mean(number+1,:) = PCA_pc([Xc(number+1);Yc(number+1)],H02,rank,U);
    pca_mean(number+2,:) = PCA_pc([Xc(number+2);Yc(number+2)],H_data(:,:,number+1),rank,U);
    number = number + 2;
end % b10 ~ b61
for i = 1:8
    pca_mean(number+1,:) = PCA_pc([Xc(number+1);Yc(number+1)],H02,rank,U);
    pca_mean(number+2,:) = PCA_pc([Xc(number+2);Yc(number+2)],H_data(:,:,number+1),rank,U);
    number = number + 2;
end % bb10 ~ bb81
for i = 1:9
    pca_mean(number+1,:) = PCA_pc([Xc(number+1);Yc(number+1)],H02,rank,U);
    pca_mean(number+2,:) = PCA_pc([Xc(number+2);Yc(number+2)],H_data(:,:,number+1),rank,U);
    number = number + 2;
end % c10 ~ c91
for i = 1:6
    pca_mean(number+1,:) = PCA_pc([Xc(number+1);Yc(number+1)],H02,rank,U);
    pca_mean(number+2,:) = PCA_pc([Xc(number+2);Yc(number+2)],H_data(:,:,number+1),rank,U);
    number = number + 2;
end % l10 ~ l61
for i = 1:6
    pca_mean(number+1,:) = PCA_pc([Xc(number+1);Yc(number+1)],H02,rank,U);
    pca_mean(number+2,:) = PCA_pc([Xc(number+2);Yc(number+2)],H_data(:,:,number+1),rank,U);
    number = number + 2;
end % r10 ~ r61
for i = 1:6
    pca_mean(number+1,:) = PCA_pc([Xc(number+1);Yc(number+1)],H02,rank,U);
    pca_mean(number+2,:) = PCA_pc([Xc(number+2);Yc(number+2)],H_data(:,:,number+1),rank,U);
    number = number + 2;
end % t10 ~ t61
for i = 1:8
    pca_mean(number+1,:) = PCA_pc([Xc(number+1);Yc(number+1)],H02,rank,U);
    pca_mean(number+2,:) = PCA_pc([Xc(number+2);Yc(number+2)],H_data(:,:,number+1),rank,U);
    number = number + 2;
end % tt10 ~ tt81
for i = 1:10
    pca_mean(number+1,:) = PCA_pc([Xc(number+1);Yc(number+1)],H02,rank,U);
    pca_mean(number+2,:) = PCA_pc([Xc(number+2);Yc(number+2)],H_data(:,:,number+1),rank,U);
    pca_mean(number+3,:) = PCA_pc([Xc(number+3);Yc(number+3)],H_data(:,:,number+2),rank,U);
    number = number + 3;
end % 11 ~ 103

save(save_doc + date + "-pca_mean","pca_mean") 
% save the pca_mean for Local_dataset_GP_evaluation.m

%% 3.1 The transition update calculation
number = 0;
for i = 1:97
    number = number + 1;
    H_data(:,:,i) = H_data(:,:,i) - H0;
end % docName1 & docName2

number = 97;
for i = 1:6
    H_data(:,:,number + 2) = H_data(:,:,number + 2) - H_data(:,:,number + 1);
    H_data(:,:,number + 1) = H_data(:,:,number + 1) - H02;
    number = number + 2;
end % b10 ~ b61
for i = 1:8
    H_data(:,:,number + 2) = H_data(:,:,number + 2) - H_data(:,:,number + 1);
    H_data(:,:,number + 1) = H_data(:,:,number + 1) - H02;
    number = number + 2;
end % bb10 ~ bb81
for i = 1:9
    H_data(:,:,number + 2) = H_data(:,:,number + 2) - H_data(:,:,number + 1);
    H_data(:,:,number + 1) = H_data(:,:,number + 1) - H02;
    number = number + 2;
end % c10 ~ c91
for i = 1:6
    H_data(:,:,number + 2) = H_data(:,:,number + 2) - H_data(:,:,number + 1);
    H_data(:,:,number + 1) = H_data(:,:,number + 1) - H02;
    number = number + 2;
end % l10 ~ l61
for i = 1:6
    H_data(:,:,number + 2) = H_data(:,:,number + 2) - H_data(:,:,number + 1);
    H_data(:,:,number + 1) = H_data(:,:,number + 1) - H02;
    number = number + 2;
end % r10 ~ r61
for i = 1:6
    H_data(:,:,number + 2) = H_data(:,:,number + 2) - H_data(:,:,number + 1);
    H_data(:,:,number + 1) = H_data(:,:,number + 1) - H02;
    number = number + 2;
end % t10 ~ t61
for i = 1:8
    H_data(:,:,number + 1) = H_data(:,:,number + 1) - H02;
    H_data(:,:,number + 2) = H_data(:,:,number + 2) - H02 - H_data(:,:,number + 1);
    number = number + 2;
end % tt10 ~ tt81
for i = 1:10
    H_data(:,:,number + 3) = H_data(:,:,number + 3) - H_data(:,:,number + 2);
    H_data(:,:,number + 2) = H_data(:,:,number + 2) - H_data(:,:,number + 1);
    H_data(:,:,number + 1) = H_data(:,:,number + 1) - H02;
    number = number + 3;
end % 11 ~ 103
% docName3

%% 3.2 The error model

% the center, volume -> error data
dep_center = [Xc,Yc];
Vol = 7.5889e+04; % approx volume of full amount
H_error = zeros(m,n,number);

if nominal_model_flag == 0
    % for function_input_2d
    lambdaX = 2.302227968708186e+03; lambdaY = 1.631169900899356e+03;
    kV = 3.445999018892074; Xmove = 2.1830e-07; Ymove = 0.207118001857057;
    Sigma = [lambdaX,0;0,lambdaY];
    for i = 1:number
        delta_H = function_input_2d(X,Y,dep_center(i,:)-[Xmove,Ymove],kV*Vol_data(i)*Vol,Sigma,the,xf,yr,yl);
        H_error(:,:,i) = H_data(:,:,i) - delta_H;   % error data
        %     close all % test use (study about the error distribution!)
        %     figure; mesh(delta_H); title("Gaussian pdf")
        %     figure; mesh(H_error(:,:,i)); title("Error Model")
    end
    
elseif nominal_model_flag == 1
    % for unit_sym_input_2d
    lambdaX = 1.209827074460249e+03; lambdaY = 1.857790332589524e+03;
    kV = 3.9410; Xmove = 0.0253; Ymove = 0.1487; alpha = 0.6147;
    Sigma = [lambdaX,0;0,lambdaY];
    for i = 1:number
        delta_H = unit_sym_input_2d(X,Y,dep_center(i,:)-[Xmove,Ymove],kV*Vol_data(i)*Vol,Sigma,alpha);
        H_error(:,:,i) = H_data(:,:,i) - delta_H;   % error data
        %     close all % test use (study about the error distribution!)
        %     figure; mesh(delta_H); title("Gaussian pdf")
        %     figure; mesh(H_error(:,:,i)); title("Error Model")
    end
end

%% 4.1 Determine the local area for each data
% Output:
% H_local(m,n,number) - the data from H_error, outside local area is all 0.
% X_local(m,n,number) - the data of X, outside local area is all 0.
% Y_local(m,n,number) - the data of Y, outside local area is all 0.

% Parameters
ux = 70; uy = 70;   % area upper bound to center
lx = 70; ly = 70;   % area lower bound to center

% Get the area of local space
X_range = [Xc-lx, Xc+ux]; Y_range = [Yc-ly, Yc+uy];
Local_area = ones(m,n,number); % 1 and 0 meshgrid, 1 represents the local area
H_local = zeros(m,n,number); % The data of local area
X_local = zeros(m,n,number);
Y_local = zeros(m,n,number); % the meshgrid of position value for local area

for k = 1:number
    % represent with logical value
    X_area = X(1,:)<X_range(k,1); X_area = X_area | X(1,:)>X_range(k,2);
    Y_area = Y(:,1)<Y_range(k,1); Y_area = Y_area | Y(:,1)>Y_range(k,2);
    
    % convert logical value to double
    for i = 1:length(Y_area)
        for j = 1:length(X_area)
            if X_area(j) == 0 && Y_area(i) == 0
                Local_area(i,j,k) = 1;
            else
                Local_area(i,j,k) = 0;
            end
        end
    end
    
    % extract for X and Y
    X_local(:,:,k) = Local_area(:,:,k) .* X; 
    Y_local(:,:,k) = Local_area(:,:,k) .* Y;
    
    % extract the local data from data, not local => 0
    H_local(:,:,k) = Local_area(:,:,k) .* H_error(:,:,k);
end

%% 5. Make it sparse!
% Because different data has different size of local area,
% here let us make it sparse first.

step_length = 13; % level of sparseness

H_local_sparse = []; X_local_sparse = []; Y_local_sparse = [];

ix = 1; jy = 1;
for i = 1:step_length:m
    for j = 1:step_length:n
        H_local_sparse(ix,jy,:) = H_local(i,j,:);
        X_local_sparse(ix,jy,:) = X_local(i,j,:);
        Y_local_sparse(ix,jy,:) = Y_local(i,j,:);
        jy = jy + 1;
    end
    ix = ix + 1;
    jy = 1;
end

%% 6. Plot to check the local area and data
% for num = 8  % the data number: 1~61 62~97 98~225
for num = [8 15 50 123 221 222]

% Evaluation: plot H_local to see whether your local area covers well!!
figure; mesh(H_error(:,:,num)); figure; mesh(H_local(:,:,num));

% Evaluation idea: Show the local area
figure; hold on;
title("Data No." + num + ", the dataset point(red) and selected point(black)")

% full data points
plot(X,Y,'LineStyle','none','Marker','.','Color',[218,165,32]/255)

% soil loading center
plot(Xc(num),Yc(num),'LineStyle','none','Marker','x','Color','[0, 0, 0]','LineWidth',1.7)

% local area bounds
plot([X_range(num,1),X_range(num,2)],[Y_range(num,1),Y_range(num,1)],'Color','[0, 0, 0]')
plot([X_range(num,1),X_range(num,2)],[Y_range(num,2),Y_range(num,2)],'Color','[0, 0, 0]')
plot([X_range(num,1),X_range(num,1)],[Y_range(num,1),Y_range(num,2)],'Color','[0, 0, 0]')
plot([X_range(num,2),X_range(num,2)],[Y_range(num,1),Y_range(num,2)],'Color','[0, 0, 0]')

% selected data points
plot(X_local_sparse(:,:,num), Y_local_sparse(:,:,num),'LineStyle','none','Marker','x','Color','[0, 0, 0]')
hold off; % (0,0) are not aovided.

end

%% 7. Make a dataset
% pick up all the selected data points in the determined local area!
% and also their normalized relative positions.
% X_data = [Xc, Yc, Xre, Yre]
% Y_data = H

[m,n,number] = size(H_local_sparse);
X_data = [];
Y_data = [];

for k = 1:number
    % get rid of the 0 elements in H_local_sparse, X_local_sparse, Y_local_sparse
    XX = X_local_sparse(:,:,k);
    YY = Y_local_sparse(:,:,k);
    HH = H_local_sparse(:,:,k);
    for i = m:-1:1
        if sum(XX(i,:))==0 && sum(YY(i,:))==0 && sum(HH(i,:))==0
            XX(i,:) = []; YY(i,:) = []; HH(i,:) = [];
        end
    end
    for j = n:-1:1
        if sum(XX(:,j))==0 && sum(YY(:,j))==0 && sum(HH(:,j))==0
            XX(:,j) = []; YY(:,j) = []; HH(:,j) = [];
        end
    end
    
    % the relative positions and normalize
    Xre = (XX - Xc(k))/170; Yre = (YY - Yc(k))/80;
    Xre = Xre(:); Yre = Yre(:);
    Xc_nor = Xc(k)/170 * ones(length(Xre),1);
    Yc_nor = Yc(k)/80 * ones(length(Yre),1);
    
    % Volume
    Vol_nor = Vol_data(k) * ones(length(Xre),1);
    
    % the PCs & mean value of loading area
    pca_mean_nor = [pca_mean(k,1:rank)/800, pca_mean(k,rank+1)/40];
    % PCs: 758.0602, mean: 37.6060
    pca_mean_nor = ones(length(Xre),1) * pca_mean_nor;
    
    % save the data to dataset vector
    Y_data = [Y_data; HH(:)];
    X_data = [X_data; Xc_nor, Yc_nor, Xre, Yre, Vol_nor, pca_mean_nor];
end

disp("The input data N = " + length(X_data(:,1)) + ", D = " + length(X_data(1,:)))

% X_test
pca_mean_vec = PCA_pc([85,0],H0,rank,U);
pca_mean_vec = ones(9400,1)*[pca_mean_vec(1:rank)/800, pca_mean_vec(rank+1)/40];
X_test = [0.5*ones(9400,1), zeros(9400,1), (X(:)-85)/170, Y(:)/80,...
    1.0*ones(9400,1),pca_mean_vec];

H0 = initializeH_2d(X,xc,xr,thr,thf);

filename = date + "-local_dataset.mat";
save(save_doc + filename,'X_data','Y_data', 'X', 'Y', 'H0')
disp("The data is saved as "+ filename)