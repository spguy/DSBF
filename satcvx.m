% basic parameters
sat_num = 4;
H = 0.1;%地面范围
step = H/201;
y = -H/2:step:H/2;
x = -H/2:step:H/2;
node_num = 10;
%node_x = randi(length(x_ground),1,node_num);
%node_y = randi(length(y_ground),1,node_num);
%好靶子
%node_x = [158	73	97	129	4	42	99	147	45	181];
%node_y = [ 195	28	163	40	43	11	17	77	20	194];
%对比组
node_x = [117	70	147	151	184	156	110	166	173	170];
node_y = [ 163	35	59	104	201	184	115	9	113	36];
d = 1.92;%satellite distance
D = 550;%distance to the ground
nvars = sat_num*2+sat_num*2;%变量个数
%satellite coordinates 
c = 3e8;
frequency = 3.5e9;
lambda = c/frequency;%波长
lambda = lambda/1000;%换算成km
%satellite coordinates 
t = 0:1e-11:1e-9;
s = zeros(sat_num,length(t));

baselinex = [zeros(1,sat_num) ones(1,sat_num) d d -d -d d -d d -d];%[phi A sat_X sat_Y]
[baseline_power] = fitness_power_nodeonly_phiAloc(baselinex,sat_num,x_ground,y_ground,node_x,node_y,D);%以无优化的平均功率作为门限
threshold = baseline_power;
xx = zeros(1,nvars);

%disp(['Debug: x = ', num2str(x)]);
phi = xx(1:sat_num);
A = xx((sat_num+1):sat_num*2);
sat_loc_x = xx(sat_num*2+1:sat_num*3);     
sat_loc_y = xx(sat_num*3+1:sat_num*4);
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

%threshold = 4; 
score_power = 0;
score_propability = 0;
for m = 1:length(node_x)
    score_power = score_power + cor(m);
    if cor(m)>abs(threshold)
        score_propability = score_propability + 1;
    end
end
score_power = -score_power/length(node_x);%求平均
power = score_power;
score_power = score_power/sat_num.^2;%归一化
score_propability = -score_propability/length(node_x);
propability = score_propability;
weight_power = 0.5;
weight_propability = 1-weight_power;
score = weight_power*score_power + weight_propability*score_propability;
cvx_begin
variable xx(nvars) 
minimize(score)
subject to 
sum(xx(sat_num+1:2*sat_num).^2) - sat_num <= 0.1;
cvx_end