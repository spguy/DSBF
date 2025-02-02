function powerdistribution(cor,x_plain,y_plain,node_x,node_y)
[X,Y] = meshgrid(x_plain,y_plain);
figure
meshc(X,Y,cor')
title('power distribution/db'); 
xlabel('x-axis','fontsize',64);  
ylabel('y-axis','fontsize',64); 
zlabel('power /W','fontsize',64)

figure
imagesc(X(1, :), Y(:, 1), cor');
hold on
plot(x_plain(node_x), y_plain(node_y), 'ro ','MarkerSize', 20, 'MarkerFaceColor', 'r');
set(gca,'FontSize',64,'FontName','Times New Roman')
legend('Ground Targets','FontSize',64)
h = colorbar;
ylabel(h,'Power Gain(dB)','fontsize',64,'FontName','Times New Roman');
xlabel('x-axis(m)','fontsize',64,'FontName','Times New Roman');  
ylabel('y-axis(m)','fontsize',64,'FontName','Times New Roman'); 
