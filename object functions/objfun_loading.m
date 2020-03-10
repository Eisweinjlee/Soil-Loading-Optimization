function J = objfun_loading(u,H0,R,Nominal_model,X_data,Y_data,X,Y,hyp_sparseGP,U)

times = length(u)/3;
H_last = H0;

% Loading soil to the end
for i = 1:times
    H_after = gp_predict(H_last,u(3*i-2),u(3*i-1),u(3*i),U,X,Y,...
        X_data,Y_data,hyp_sparseGP,Nominal_model);
    H_last = H_after;
end

% evaluation the result
J = immse(H_after,R);

end