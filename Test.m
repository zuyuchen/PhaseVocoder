syms a a1 a2 a3 a4
P = ((a - a2)*(a-a3)*(a-a4))/((a1-a2)*(a1-a3)*(a1-a4));
[C, T] = coeffs(expand(P), a);
