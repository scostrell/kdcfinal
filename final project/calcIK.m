function [iter,thetas] = calcIK(twists, theta_init, g_init, pos)

    pos_des = pos(1:3)';
    quat_des = quaternion(pos(end),pos(4),pos(5),pos(6));
    theta_curr = theta_init;
    err_norm = 1;
    iter = 0;
%     kp = 1;
%     ko = 1;
    kp = 0.25;
    ko = 1;
    timestep = 0.002;
    thetas = [theta_init];

    while((err_norm > 1e-3) && (iter < 5000))
        iter = iter + 1;
        g = calcG(twists,theta_curr,g_init);
        pos_curr = g(1:3,4);
        quat_curr = quaternion(g(1:3,1:3),'rotmat','point');

        pos_err = pos_des - pos_curr;
        quat_err = quat_des*(quat_curr.conj);
        quat_err_compact = quat_err.compact';

        err_norm = norm([pos_err; quat_err_compact(1)*quat_err_compact(2:end)]);

        V = [kp*pos_err; ko*quat_err_compact(1)*quat_err_compact(2:end)];
        J = calcJacobian(twists,theta_curr);
        theta_dot = pinv(J)*V;
        theta_curr = theta_curr + theta_dot*timestep;
        thetas = [thetas,theta_curr];
    end
%     thetas = theta_curr;
end