%% Record the coordinates of adjacent tissues
function side = str2nor_tum_side(m, n, x, y, border)

    side = zeros(1, 4);
    if m - border > 0
        side(1) = m - border;
    else
        side(1) = 1;
    end
    if m + border < x
        side(2) = m + border;
    else
        side(2) = x;
    end
    if n - border > 0
        side(3) = n - border;
    else
        side(3) = 1;
    end
    if n + border < y
        side(4) = n + border;
    else
        side(4) = y;
    end

end