function ke = el1_ke(xe, properties)

Le = abs(xe(2) - xe(1));
E = properties.E;
A = properties.A;

ke = E * A / Le * [+1.0, -1.0;
    -1.0, +1.0];

end