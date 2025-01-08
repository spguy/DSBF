clc
clear
close all
figure
for k = 1:10
track_theta = [0 90];
track_bias = [0 0];


%初始坐标为（0，1）
vs = 7.8;
update = 0.05;
sat_x = update*k*vs*cosd(track_theta);
sat_y = track_bias+update*k*vs*sind(track_theta);



frequency = 3.5e9;
c = 3e8;
lambda = c/frequency;
lambda = lambda/1000;%km

h = 550; 

ground_x = -1:0.01:1;
ground_y = -1:0.01:1;
[X,Y] = meshgrid(ground_x,ground_y);
phiplot = zeros(length(ground_y));
for i =1:length(sat_y)
distance = (sat_x(i) - X).^2 + (sat_y(i) - Y).^2 + h^2;
distance = sqrt(distance);
phiplot = phiplot + 2*pi*distance/lambda;
phiplot = mod(phiplot, 2*pi);
end



meshc(X,Y,phiplot')
title('power distribution/db'); 
xlabel('x-axis','fontsize',14);  
ylabel('y-axis','fontsize',14); 
zlabel('power /W','fontsize',14)
% Adjust view for better visualization
    view(0, 90)
    
    % Pause for a short time to create a dynamic effect
    pause(0.5)
    
    % Use drawnow to update the figure
    drawnow
    
    % Clear the current plot for the next iteration
    clf
end
