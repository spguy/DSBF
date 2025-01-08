function [c,ceq] = constraint_phiAloc(x,sat_num)
ceq = sum(x(sat_num+1:2*sat_num).^2) - sat_num;
sat_loc_x = xx(sat_num*2+1:sat_num*3);     
sat_loc_y = xx(sat_num*3+1:sat_num*4);
c = [];