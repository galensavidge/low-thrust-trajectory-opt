function r = position_from_MEE(x)

r = zeros(3, length(x));
for i = 1:length(r)
    % Translate MEE to classical orbital elements (COE)
    coe = MEE2COE(x(:,i));

    % Extract orbital elements from state vectors
    p = x(1,i);
    e = coe(2);
    inc = coe(3);
    RAAN = coe(4);
    AOP = coe(5);
    true_anom = coe(6);
    
    % Construct position in the local orbital (LVLH) frame
    r_norm = p/(1 + e*cosd(true_anom));
    r_LVLH = [r_norm; 0; 0];
    
    % Express in inertial frame
    C = C_LVLH_to_inertial(RAAN, inc, AOP, true_anom);
    r(:,i) = C*r_LVLH;
end