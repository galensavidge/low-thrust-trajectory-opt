function x_soln = newton_iteration(f, x0, scaling_factor, delta, ...
    max_iter, min_step)

min_step_2 = min_step^2;
num_iter = 0;
x = x0;
while true
    G = finite_difference_jacobian(f, x, delta);
    step = scaling_factor*inv(G)*f(x);
    x = x - step
    
    if step.^2 < min_step_2
        break
    end

    num_iter = num_iter + 1;
    if num_iter > max_iter
        warning('Newton iteration exiting because num_iter > max iter')
        break
    end
end

x_soln = x;