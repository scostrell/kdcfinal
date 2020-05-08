function jac = calcJacobian(twists,thetas)
% a function using twists to find the jacobian (see MLS pp. 116-117)
    jac = zeros(6,numel(thetas));
    jac(:,1) = twists(:,1);
    for i = 2:numel(thetas)
        exps = eye(4); 
        for j = 1:i
            exps = exps*twistExp(twists(:,j),thetas(j));
        end
        jac(:,i) = adjoint(exps)*twists(:,i);
    end
end