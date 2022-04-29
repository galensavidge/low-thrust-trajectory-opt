function [best_p0, xf_err] = solve_2pbvp(P, t0, tf, rho, p0_guess)
% Solves the two-point boundary value problem posed by indirect single
% shooting, in this case using the modified equinoctial elements.

shot_wrapper = @(p0)(single_shot(P, t0, tf, p0(1:6), p0(7), rho));

options = optimset('FinDiffRelStep', 1e-6, 'Display', 'iter-detailed', ...
    'TolX', 1e-15, 'UseParallel', true);
[best_p0, xf_err] = fsolve(shot_wrapper, p0_guess, options);

end

function xf_constraint = single_shot(P, t0, tf, px0, pm0, rho)
    [~, X_hist] = propagator_MEE_indirect(P, t0, tf, px0, pm0, rho);
    Xf = X_hist(:,end);
    
    % Constrain final orbit state = xf, pL_f = 0, and pm_f = -1
    xf_constraint = [Xf(1:5) - P.xf(1:5); Xf(13); Xf(14) + 1];
end