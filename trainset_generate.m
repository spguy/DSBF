% basic parameters
sat_num = 4;
node_loc = zeros(node_num*2,1);
phiAloc = zeros(sat_num*4,1);
H = 1;%地面范围
step = H/2001;
y_ground = -H/2:step:H/2;
x_ground = -H/2:step:H/2;    
D = 550;%distance to the ground
nvars = sat_num*2+sat_num*2;%变量个数
d = 5; 
baselinex = [zeros(1,sat_num) ones(1,sat_num) d d -d -d d -d d -d];%[phi A sat_X sat_Y]
tic
for i = 1:50000
disp("训练集个数")
disp(i)
node_x = randi(length(x_ground),1,node_num);
node_y = randi(length(y_ground),1,node_num);
nodes = [x_ground(node_x),y_ground(node_y)]';
% 固定卫星功率
[baseline_power] = fitness_power_nodeonly_phiAloc(baselinex,sat_num,x_ground,y_ground,node_x,node_y,D);%以无优化的平均功率作为门限
threshold = baseline_power;
targetFunction = @(x)fitness_weighted_nodeonly_phiAloc(x,sat_num,x_ground,y_ground,node_x,node_y,D,threshold); % 目标函数：调用 fitness 函数
constraintFunction = @(x)constraint(x,sat_num);
UB =[pi*ones(1,sat_num) 2*ones(1,sat_num) d 0 0 d d d 0 0];%卫星分别位于四个象限
LB =[-pi*ones(1,sat_num) zeros(1,sat_num) 0 -d -d 0 0 0 -d -d];
% 使用 ga 函数进行单目标优化
[x, fval, exitflag, output,population,scores] = ga(targetFunction, nvars, [], [], [], [], LB, UB,constraintFunction);
node_loc = cat(2,node_loc,nodes);
phiAloc = cat(2,phiAloc,x');

end
toc