function J = objfun_loading_symAtY(u,H0,R,Nominal_model,X_data,Y_data,X,Y,hyp_sparseGP,U)

loading_times = length(u)/2;
H_last = H0;

% Loading soil to the end
for i = 1:loading_times
    H_after = gp_predict(H_last,u(2*i-1),0,u(2*i),U,X,Y,...
        X_data,Y_data,hyp_sparseGP,Nominal_model);
    H_last = H_after;
end

% evaluation the result
J = immse(H_after,R);

end