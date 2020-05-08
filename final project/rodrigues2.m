function exp_vec = rodrigues2(omega,theta)
% vector exponential for non unit omega p. 28 MLS
    omega_hat = matcross(omega);
    exp_vec = eye(3) + omega_hat/norm(omega)*sin(norm(omega)*theta) + ...
        omega_hat^2/norm(omega)^2*(1-cos(norm(omega)*theta));
end