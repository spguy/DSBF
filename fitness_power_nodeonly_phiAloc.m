function [score] = fitness_power_nodeonly_phiAloc(xx,sat_num,x,y,node_x,node_y,D)
phi = xx(1:sat_num);
A = xx((sat_num+1):sat_num*2);
sat_loc_x = xx(sat_num*2+1:sat_num*3);     
sat_loc_y = xx(sat_num*3+1:sat_num*4);
%satellite coordinates 
c = 3e8;
frequency = 3.5e9;
lambda = c/frequency;%波长
lambda = lambda/1000;%换算成km
%satellite coordinates 
t = 0:1e-11:1e-9;
s = zeros(sat_num,length(t));
for n = 1:sat_num
    s(n,:) = A(n)*exp(1i*(2*pi*frequency*t+phi(n)));
end
cor = zeros(1,length(node_x));
R = zeros(1,sat_num);
delta_phi = zeros(1,sat_num);
signal_receive = zeros(1,length(t));
%%
   for m = 1:length(node_x)
        for k = 1:sat_num
            R(k) = sqrt((y(node_y(m))-sat_loc_y(k)).^2+(x(node_x(m))-sat_loc_x(k)).^2+D^2);
            delta_phi(k) = R(k)*2*pi/lambda;
            signal_receive = signal_receive + s(k,:)*exp(1i*delta_phi(k));
        end
        cor(m) = signal_receive*signal_receive'/length(signal_receive);
        signal_receive = zeros(1,length(t));
   end
score = 0;
for m = 1:length(node_x)
        score = score + cor(m);
end
score = -score/length(node_x);
end