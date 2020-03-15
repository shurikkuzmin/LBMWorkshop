%% Contour plotting
clc
clear all

phase=dlmread('phase020000.dat',' ');
%imagesc(phase)
figure(1)
[c,h]=contour(phase,1)
hold on

aver=mean(c,2)

c=c(:,2:end)
crel=sqrt((c(1,:)-aver(1)).^2+(c(2,:)-aver(2)).^2)
rmin=min(crel)
rmax=max(crel)

x0=aver(1);
y0=aver(2);

xaxis=find(phase(round(x0),:)>0)
yaxis=find(phase(:,round(y0))>0)

radx=max(yaxis)-y0;
rady=max(xaxis)-x0;
theta=-atan(rady/radx)

a=rmin;
b=rmax;
alpha=0:pi/180:2*pi;
x=a*cos(alpha);
y=b*sin(alpha);

xprime=x0+x.*cos(theta)+y.*sin(theta);
yprime=y0-x.*sin(theta)+y.*cos(theta);
%plot(xprime,yprime,'r+')

d=abs(rmax-rmin)/(rmax+rmin)
%% Postprocessing with the Matlab


%plot(phase(dims(1)/2,:),'r+')
%hold on
%plot(phase(:,(dims(2)-1)/2),'bo')