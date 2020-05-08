function C = coriolis(M, th1, th2, th3, dth1, dth2, dth3)
% Using eqn 4.23 on p.170 of MLS, assuming robot in R3
    C = sym('C',[3 3]); % make symbolic matrix to avoid sym vs double error
    Csym = C;
    theta = [th1, th2, th3];
    dtheta = [dth1, dth2, dth3];
    for i = 1:3
        for j = 1:3
            for k = 1:3 
                % see diff docs
                C(i,j) = C(i,j) + (diff(M(i,j), theta(k)) + diff(M(i,k), theta(j))...
                    - diff(M(k,j),theta(i)))*(0.5*dtheta(k));
            end
        end
    end
    C = C - Csym;
end