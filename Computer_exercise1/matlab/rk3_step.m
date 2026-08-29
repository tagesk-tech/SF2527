function unext = rk3_step(f, tn, un, h)
    k1 = f(tn, un);
    k2 = f(tn + h, un + h*k1);
    k3 = f(tn + h/2, un + h*k1/4 + h*k2/4);

    unext = un + h*(k1 + k2 + 4*k3)/6;
end
