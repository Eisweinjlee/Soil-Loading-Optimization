%% make prediction
% Input: H0, xc, yc, Vol
% Output: H

function H_after = gp_predict(H_last,Xc,Yc,Vol,U,X,Y,X_data,Y_data,hyp_sparseGP,Nom)
rank = 4;
Xc_nor = Xc/170; Yc_nor = Yc/80;

% pca_mean normalization
pca_mean_vec = PCA_pc([Xc,Yc],H_last,rank,U);
pca_mean_vec = ones(9400,1)*[pca_mean_vec(1:4)/800, pca_mean_vec(4+1)/40];

% query data x*
X_ast = [Xc_nor*ones(9400,1), Yc_nor*ones(9400,1), ...
    (X(:)-Xc)/170, (Y(:)-Yc)/80, Vol*ones(9400,1), pca_mean_vec];

% configuration of GP prediction
xu = hyp_sparseGP.xu; cov = {'apxSparse', {@covSEard}, xu};
infmethod = @infGaussLik;
inff = @(varargin) infmethod(varargin{:},struct('s', 0)); 

% prediction
[ymu,ys2] = gp(hyp_sparseGP, inff, [], cov, {@likGauss},...
    X_data, Y_data, X_ast);

% reshape the result
m = 94; n = 100;
H_error_pred = reshape(ymu,[m,n]);  % mean of GPs
CovDist = reshape(ys2,[m,n]);       % covariance of GPs

%% 3. Relative covariance modification

% find the data close to center
disX = (X(1,:)-Xc).^2; disY = (Y(:,1)-Yc)'.^2;
[~,minX_p] = min(disX); [~,minY_p] =  min(disY);
% data_cen = [X(1,minX_p);Y(minY_p,1)];

% calculate all relative covariance
cov_cen = CovDist(minY_p,minX_p); % the covariance of center
w = cov_cen ./ CovDist;

% Sigmoid function
a = 1; b = 18.1; % https://www.desmos.com/calculator/kn9tpwdan5
gw = 1./(1 + exp(a + -b.*w)); % sigmoid

% modify the mean predition through sigmoid
H_modified = gw .* H_error_pred;

%% 4. Add to the nominal Gaussian pdf model
excavator_data;
bucket_vol = 7.5889e+04;

if Nom(1) == 0
    % Nominal_model = [nominal_model_flag,Xmove,Ymove,kV,Sigma,the,xf,yr,yl];
    H_nominal = function_input_2d(X,Y,[Xc;Yc]-[Nom(2);Nom(3)],Nom(4)*Vol*bucket_vol,[Nom(5),0;0,Nom(6)],Nom(7),Nom(8),Nom(9),Nom(10));
elseif Nom(1) == 1
    % Nominal_model = [nominal_model_flag,Xmove,Ymove,kV,Sigma;alpha];
    H_nominal = unit_sym_input_2d(X,Y,[Xc;Yc]-[Nom(2),Nom(3)],Nom(4)*Vol*bucket_vol,[Nom(5),0;0,Nom(6)],Nom(7));
end

H_noGP = H_last + H_nominal;
H_after = H_noGP + H_modified;

end
