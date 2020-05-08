function C = coriolis_barett(M, th1, th2, th3, th4, th5, th6, th7, dth1, dth2, dth3, dth4, dth5, dth6, dth7)
% Using eqn 4.23 on p.170 of MLS, assuming robot in R3
    n = 7;
    C = sym(zeros(n,n)); % symbolic zeros https://www.mathworks.com/matlabcentral/answers/162423-how-to-create-on-an-efficient-way-zeros-in-a-symbolic-matrix
    theta = [th1, th2, th3, th4, th5, th6, th7];
    dtheta = [dth1, dth2, dth3, dth4, dth5, dth6, dth7];
    for i = 1:n
        for j = 1:n
            for k = 1:n 
                % see diff docs
                C(i,j) = C(i,j) + (diff(M(i,j), theta(k)) + diff(M(i,k), theta(j))...
                    - diff(M(k,j),theta(i)))*(0.5*dtheta(k));
            end
        end
    end
    
end