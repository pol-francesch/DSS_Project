function R = get_hill_frame(rv)
    r = rv(1:3);
    v = rv(4:6);

    h = cross(r,v);

    r_hat = r / norm(r);
    h_hat = h / norm(h);
    y_hat = cross(h_hat, r_hat);

    R = [r_hat'; y_hat'; h_hat'];
end

