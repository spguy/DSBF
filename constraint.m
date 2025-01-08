function [c,ceq] = constraint(x,sat_num)
ceq = sum(x(sat_num+1:2*sat_num).^2) - sat_num;
c = [];