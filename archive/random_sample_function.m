function best_x = random_sample_function(f, xmin, xmax, num_samples)
% Attempts to find an approximate solution for vector-valued function f 
% between [xmin, xmax] by random sampling

N = length(xmin);

best_f_norm = Inf;

disp('Starting random sampling...')

for i=1:num_samples
    x = xmin + rand(N, 1).*(xmax - xmin);
    norm_f_2 = sum(f(x).^2);
    if norm_f_2 < best_f_norm
        best_x = x;
        best_f_norm = norm_f_2;
    end

    if rem(i, 25) == 0
        disp(string(i))
    end
end