function J = objfun_3times_search(u,H0,R,Nominal_model,X_data,Y_data,X,Y,hyp_sparseGP,U)

H_last = H0;

for i = 1:3
    if u(3*i-2) < 43
        u(3*i-2) = 43;
    elseif u(3*i-2) > 127
        u(3*i-2) = 127;
    end
    
    if u(3*i-1) < -39
        u(3*i-1) = -39;
    elseif u(3*i-1) > 39
        u(3*i-1) = 39;
    end
    
    if u(3*i) < 0
        u(3*i) = 0;
    elseif u(3*i) > 1
        u(3*i) = 1;
    end
end

% Loading soil to the end
for i = 1:3
    H_after = gp_predict(H_last,u(3*i-2),u(3*i-1),u(3*i),U,X,Y,...
        X_data,Y_data,hyp_sparseGP,Nominal_model);
    H_last = H_after;
end

% evaluation the result
J = immse(H_after,R);

end