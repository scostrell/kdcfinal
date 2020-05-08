function exp_vec = rodrigues(omega,theta)
% vector exponential, p. 28 MLS
    omega = omega/norm(omega);
    omega_hat = matcross(omega);
    exp_vec = eye(3) + omega_hat*sin(theta) + omega_hat^2*(1-cos(theta));
end