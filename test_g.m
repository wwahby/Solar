function out = test_g(W, m_vec, Npts_x)

colors = {'k', 'b', 'r', 'g', 'm'};

cind = 0;
figure(1)
clf
hold on
for m_max = m_vec
    cind = cind + 1;
    col_ind = mod(cind-1, length(colors)) + 1;
    [g, x] = calc_g(W, m_max, Npts_x);
    plot(x, g, 'color', colors{col_ind})
end

out = 0;
end
    

function [g, x]= calc_g( W, m_max, Npts_x)

x = linspace(-W, W, Npts_x);

g = accumulate_gvecs(W, m_max, x, Npts_x);


end





function gacc = accumulate_gvecs( W, m_max, x, Npts_x )
    gvec = @(m) -4*(-1)^m/(2*m-1)/pi * cos( (2*m-1)* pi*x/W);
    gacc = zeros(1,Npts_x);
    for m = 1:m_max
        gacc = gacc + gvec(m);
    end
end