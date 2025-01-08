%随机生成节点，测试平均功率和概率
function [average_power,average_propability] = node_test(xx,sat_num,sat_loc_x,sat_loc_y,x,y,node_num,D,threshold)
test_time = 100;
average_power_rec = [];
average_propability_rec = [];
for i = 1:test_time
    node_x = randi(length(x),1,node_num);
    node_y = randi(length(y),1,node_num);
    phi = xx(1:sat_num);
    A = xx((sat_num+1):sat_num*2);
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
    score_power = 0;
    score_propability = 0;
    for m = 1:length(node_x)
        score_power = score_power + cor(m);
        if cor(m)>abs(threshold)
            score_propability = score_propability + 1;
        end
    end
    average_power_rec = [average_power_rec -score_power/length(node_x)];
    average_propability_rec = [average_propability_rec -score_propability/length(node_x)];
end
average_power = -mean(average_power_rec);
average_propability = -mean(average_propability_rec);
end