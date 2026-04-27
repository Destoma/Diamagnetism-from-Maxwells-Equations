function integral = comp_simpson(func_handle, num_intervals, x_start, x_end)
    f = func_handle;
    endpts = linspace(x_start, x_end, num_intervals + 1);
    interval_length = (x_end - x_start) / num_intervals;
    simp_values = zeros(num_intervals);

    for i = 1:num_intervals
        simp_values(i) = (interval_length / 6) * (f(endpts(i)) + 4*f(endpts(i) + (interval_length / 2)) + f(endpts(i+1)));
    end

    integral = 0;
    for i = 1:num_intervals
        integral = integral + simp_values(i);
    end
end