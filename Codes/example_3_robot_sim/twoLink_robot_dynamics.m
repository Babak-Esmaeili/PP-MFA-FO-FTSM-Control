function qddot = twoLink_robot_dynamics(t,q,qdot,u,d)

    g = 9.81;   % m/(s^2)

    L1 = 1;   % m
    L2 = 1;   % m
%     L1 = 0.5;   % m
%     L2 = 0.5;   % m
    
    m1 = 1 + 0.1*sin(t);     % kg
    m2 = 1 + 0.1*sin(t);     % kg

    J1 = 0.1;     % kg/(m^2)
    J2 = 0.1;     % kg/(m^2)

    m_11 = m1*L1 + m2*(L1^2 + L2^2 + 2*L1*L2*cos(q(2))) + J1 + J2;
    m_12 = m2*(L2^2 + L1*L2*cos(q(2))) + J2;
    m_21 = m_12;
    m_22 = m2*L2^2 + J2;

    M = [ m_11  m_12
          m_21  m_22 ];

    c_11 = -m2*L1*L2*qdot(2)*sin(q(2));
    c_12 = -m2*L1*L2*(qdot(1)+qdot(2))*sin(q(2));
    c_21 = m2*L1*L2*qdot(1)*sin(q(2));
    c_22 = 0;
     
    C = [ c_11  c_12
          c_21  c_22 ];

    gamma1 = (m1+m2)*L1*cos(q(1)) + m2*L2*cos(q(1)+q(2));
    gamma2 = m2*L2*cos(q(1)+q(2));
    
    G = [ gamma1*g
          gamma2*g ];

    qddot = M^-1*(-C*qdot-G) + (M^-1)*u + (M^-1)*d;

end
