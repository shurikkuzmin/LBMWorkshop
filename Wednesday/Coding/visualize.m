%% The construction of movie from data files
clear all
clc

a=dir('height*.dat')
for counter=1:size(a)
    file_name=a(counter).name;
    rho=dlmread(file_name);
    %surf(rho)
    %zlim([0,3.0]);
    imagesc(rho)
    F(counter)=getframe;
end
movie(F)
