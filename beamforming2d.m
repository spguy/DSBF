function [gain] = beamforming2d(M,N,theta,phi,frequency,angle)
%initializin
%M = 16;
%N = 16;%array size=M*N
%theta = -10;%方位角
%phi = -10;%俯仰角
theta_all = angle;
phi_all = angle;
c = 3e8;
wavelength = c/frequency;
dx = wavelength/2;
dy = wavelength/2;
kx = 2*pi/wavelength*dx;
ky = 2*pi/wavelength*dy;

%steering vector 
    ax = exp(-1j*kx*(0:M-1)*sind(theta)*sind(phi));
    ay = exp(-1j*ky*(0:N-1)*cosd(theta)*sind(phi));
    axy = ax.'*ay;
 %%
% %directional gain
gain = zeros(length(theta_all),length(phi_all));
for i = 1:length(theta_all)
    for j = 1:length(phi_all)
        axx = exp(1j*kx*(0:M-1)*sind(theta_all(i))*sind(phi_all(j)));
        ayy = exp(1j*ky*(0:N-1)*cosd(theta_all(i))*sind(phi_all(j)));
        axxyy = axx.'*ayy;
        gain(i,j) = sum(sum(axxyy.*axy)); 
    end
end

[X,Y] = meshgrid(theta_all,phi_all);
meshc(X,Y,abs(gain))
axis([-90 90 -90 90]);
xlabel('俯仰角');
ylabel('方位角');

end