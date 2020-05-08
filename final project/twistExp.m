function exp_twist = twistExp(twist,theta)
% twist exponential, p.42 MLS
    v = twist(1:3);
    omega = twist(4:6);

    exp_twist = [rodrigues2(omega, theta), (eye(3) - rodrigues2(omega,theta))...
        *cross(omega,v) + omega*dot(omega,v)*theta; 0 0 0 1];
end