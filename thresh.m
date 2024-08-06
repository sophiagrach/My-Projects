function result = thresh(f_, a_)
%
abs_f = abs(f_);
idx_a = find(abs_f >= a_);
result = zeros(length(f_), 1);
result(idx_a) = f_(idx_a);
return
