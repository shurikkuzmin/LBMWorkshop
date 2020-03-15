%% Visualization (uncomment for nice movie)
%a=dir('phase*.dat')
%figure(1)
%for counter=1:size(a)
%    file_name=a(counter).name;
%    rho=dlmread(file_name);
%    %surf(rho)
%    imagesc(rho)
%    zlim([-1.1,1.1]);
%    F(counter)=getframe;
%end

%% Getting countour (wait until the second figure appears)
figure(2)
rho=dlmread('phase020000.dat');
rho=flipud(rho);
[c,h]=contour(rho,1);
c=c(:,2:end);

%[coormax,indmax]=max(c(1,:))
%[coormin,indmin]=min(c(1,:))

coeff=(c(2,2)-c(2,1))/(c(1,2)-c(1,1))
theta=atan(coeff)
atan(coeff)*180/pi

wall_gradient=0.15
a=0.04
k=0.04
x = fzero(@(x) sqrt(2)*sqrt(cos(acos(sin(x)^2)/3)*(1-cos(acos(sin(x)^2)/3)))-wall_gradient,theta)

param=cos(acos(sin(theta)^2)/3);
wall_gradient_computed=sqrt(2)*sqrt(param*(1-param))

Sim=theta*180/pi
Theor=x*180/pi