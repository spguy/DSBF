clc
clear
close all

NP = 50; %种群数量
%[sat1_x sat1_y sat2_x sat2_y ...]
bit = 4;
sat_num = 3;
L = bit * sat_num *2 ; %二进制位串长度
%f = round(rand(NP,L)); %随机获得初始种群
f = zeros(NP,L/2);
for sa = 1:sat_num
        f = [f ones(NP,3) zeros(NP,1)];%左为低位
end
phi = zeros(1,sat_num);
A = zeros(1,sat_num);%信号幅度
d = 1.92;%satellite distance
D = 550;%distance to the ground
sat_x = [d d -d -d];
sat_y = [d -d d -d];
%satellite coordinates 
H = 10;
step = H/201;
y_plain = -H/2:step:H/2;
x_plain = -H/2:step:H/2;
Pc = 0.8; %交叉率
Pm = 0.2; %变异率
G = 100; %最大遗传代数
PHImax = pi/2;
PHImin = -pi/2;
%总功率为4
Amax = 2;
Amin = 0;

%% 随机生成地面节点
node_num = 20;
node_x = randi(length(x_plain),1,node_num);
node_y = randi(length(y_plain),1,node_num);
%% 记录功率分布
cor_rec = zeros(G,length(x_plain),length(y_plain));
cor_temp = zeros(NP,length(x_plain),length(y_plain));
%% 遗传
phiBest_rec = [];
ABest_rec = [];
for k = 1:G
    for i = 1:NP
            U = f(i,:);
            %将个体解码为各个卫星的位置
            for kk = 1:sat_num  
                for c = 1:bit
                    phi(kk) = U(bit*(kk-1)+c)*2^(c-1) + phi(kk);
                    A(kk) = U(bit*sat_num+bit*(kk-1)+c)*2^(c-1) + A(kk);
                end
            end
            %十进制数转位置坐标
            phisig(:,i) = PHImin+(phi+1)*(PHImax-PHImin)/2^bit;
            Asig(:,i) = Amin+(A+1)*(Amax-Amin)/2^bit;
            [Fit(i),cor_temp(i,:,:)] = power_calculation(sat_num,sat_x,sat_y,phisig(:,i),Asig(:,i),x_plain,y_plain,H,node_x,node_y,D);%生成功率分布和分数
            if abs(sum(Asig(:,i).^2)-sat_num)>0.1
                Fit(i) = 0;
            end
            phi = zeros(sat_num,1);
             A = zeros(sat_num,1);
    end
    %%
    averagepower_best = max(Fit); %最大值
    averagepower_worst = min(Fit); %最小值
    rr = find(Fit==averagepower_best);
    cor_rec(k,:,:) = cor_temp(rr(1),:,:);%记录最优功率分布
    fBest = f(rr(1),:); %历代最优个体
    %sprintf('第%d次，%f,%f',k,Fit(rr(1)),averagepower_best)
    phiBest = phisig(:,rr(1));
    ABest = Asig(:,rr(1));
    phiBest_rec = [phiBest_rec phiBest];
    ABest_rec =  [ABest_rec ABest];
    if averagepower_best-averagepower_worst ~=0
        Fit = (Fit-averagepower_worst)/(averagepower_best-averagepower_worst); %归一化适应度值
    else
        Fit = ones(1,NP);
    end
    %%%%%%%%%%%%%%基于轮盘赌的选择操作%%%%%%%%%%%%%
    sum_Fit = sum(Fit);
    fitvalue = Fit./sum_Fit;%适应度
    fitvalue = cumsum(fitvalue);
    ms = sort(rand(NP,1));
    fiti = 1;
    newi = 1;
    while newi <= NP
        if (ms(newi)) < fitvalue(fiti)
            nf(newi,:) = f(fiti,:);
            newi = newi+1;
        else
            fiti = fiti+1;
        end
    end
    %%%%%%%%%%%%%%%基于概率的交叉操作%%%%%%%%%%%%%
    for i = 1:2:NP
        p = rand;
        if p < Pc
            q = rand(1,L);
            for j = 1:L
                if q(j)==1;
                    temp = nf(i+1,j);
                    nf(i+1,j) = nf(i,j);
                    nf(i,j) = temp;
                end
            end
        end
    end
    %%%%%%%%%%%%%基于概率的变异操作%%%%%%%%%%%%%%
    i = 1;
    while i <= round(NP*Pm)
        h = randi([1,NP],1,1); %随机选取一个需要变异的染色体
        for j = 1:round(L*Pm)
            g = randi([1,L],1,1); %随机选取需要变异的基因数
            nf(h,g) =~ nf(h,g);
        end
        i = i+1;
    end
    f = nf;
    f(1,:) = fBest; %保留最优个体在新种群中
    trace(k) = averagepower_best; %历代最优适应度
end
%%
figure
plot(trace)
xlabel('迭代次数')
ylabel('目标函数值')
title('平均功率进化曲线')

%%
[X Y] = meshgrid(x_plain,y_plain); 
figure
subplot(2,2,1)
meshc(X,Y,squeeze(cor_rec(1,:,:)))
title('第一次')
subplot(2,2,2)
meshc(X,Y,squeeze(cor_rec(round(G/3),:,:)))
title(sprintf('第%d次', round(G/3)));
subplot(2,2,3)
meshc(X,Y,squeeze(cor_rec(round(G/1.5),:,:)))
title(sprintf('第%d次', round(G/1.5)));
subplot(2,2,4)
meshc(X,Y,squeeze(cor_rec(G,:,:)))
title(sprintf('第%d次', G));
%%
figure
imagesc(X(1, :), Y(:, 1), squeeze(cor_rec(round(G),:,:)));
hold on
plot(x_plain(node_x), y_plain(node_y), 'ro', 'MarkerSize', 5, 'MarkerFaceColor', 'r');
title('平面图')
%%
[baseline basecor] = power_calculation(sat_num,sat_x,sat_y,[0 0 0],[1 1 1],x_plain,y_plain,H,node_x,node_y,D);
sprintf('baseline功率：%f 最优功率：%f',baseline,averagepower_best)
figure
imagesc(X(1, :), Y(:, 1), basecor);
hold on
plot(x_plain(node_x), y_plain(node_y), 'ro', 'MarkerSize', 5, 'MarkerFaceColor', 'r');
title('无优化')