function M = manip(m1, m2, m3, J1, J2, J3, Ix1, Iy1, Iz1, Ix2, Iy2, Iz2, Ix3, Iy3, Iz3)
% Calculates manipulator inertia matrix in terms of link Jacobians, via
% eqn 4.19 on p. 168 of MLS, assuming robot in R3 and principal axes chosen
% (see MLS p. 172) - outputs a 1*9 vector for computational reasons
    M1 = [eye(3)*m1, zeros(3);
        zeros(3), [Ix1, 0, 0; 0, Iy1, 0; 0, 0, Iz1]];
    M2 = [eye(3)*m2, zeros(3);
        zeros(3), [Ix2, 0, 0; 0, Iy2, 0; 0, 0, Iz2]];
    M3 = [eye(3)*m3, zeros(3);
        zeros(3), [Ix3, 0, 0; 0, Iy3, 0; 0, 0, Iz3]];
    
    M = J1'*M1*J1 + J2'*M2*J2 + J3'*M3*J3;
%     M = reshape(M',1,9);
end