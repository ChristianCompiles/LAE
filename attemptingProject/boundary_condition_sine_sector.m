function scalar = boundary_condition_sine_sector(t)
    if 1 <= t && t <= 3
        % scalar = sin(2 * pi * (t - 1)); % works fine
        % scalar = cos(1.5 * pi * (t - 1)); % discontinuous at the boundary with Gibbs fenom
        scalar = sin(2 * pi * (t - 1)) * cos(pi * (t - 1)); % works fine but a bit funkier
    else
        scalar = 0;
    end
end