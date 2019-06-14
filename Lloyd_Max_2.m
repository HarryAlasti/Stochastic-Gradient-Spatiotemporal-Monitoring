function L = Lloyd_Max(p,xx, M)

x0 = xx(1);
xend = xx(end);
y = linspace(x0, xend, M+1);
L = zeros(M, 1);
delta = xx(2) - xx(1);
for rept = 1 : 100
    for k = 1: M
        [Min, indx_inf] = min(abs(xx - y(k)));
        [Max, indx_sup] = min(abs(xx - y(k+1)));
        Sum_num = 0;
        Sum_den = 0;
        for kk = indx_inf : indx_sup
            Sum_num = Sum_num + xx(kk) * p(kk) * delta;
            Sum_den = Sum_den + p(kk) * delta;
        end,
        L(k) = Sum_num / Sum_den;
    end,
    L_left = [x0 ; L];
    L_right = [L ; xend];
    y = (L_left + L_right) / 2;
    %y(1) = xx(1);
    %y(end) = xend;
   % figure (1); hold on; plot(L);
end
return,
