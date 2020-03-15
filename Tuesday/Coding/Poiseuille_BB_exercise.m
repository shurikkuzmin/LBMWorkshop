clear all
clc
tic
% Poiseuille flow driven by constant body force

% Lattice parameters
weights=[4/9 1/9 1/9 1/9 1/9 1/36 1/36 1/36 1/36];
cx=[0 1 0 -1 0 1 -1 -1 1];
cy=[0 0 1 0 -1 1 1 -1 -1];

% Numerical parameters
NX=3;  % Number of grids points along x
NY=16;  % Number of grid points along y
NPOP=9; % Number of populations used in velocity space discretization
NSTEPS=10000;    % Number of time steps/iterations


% Simulation parameters
y_bottom=0.5; % location of bottom wall
y_top=NY+0.5; % location of top wall

Re=10;  % Reynolds number
omega=32/(20+sqrt(208));   % Relaxation frequency
kvisc=1/3*(1/omega-0.5); % Kinematic viscosity
umax=Re*kvisc/((y_top-y_bottom)) ;% umax=0.001; % Mach number (can be understood as a CFL number)


% Macroscopic parameters
rho=ones(NX,NY);
ux=zeros(NX,NY);
uy=zeros(NX,NY);

forcex=8.*umax*kvisc./((y_top-y_bottom).^2);
forcey=0;



% Initialize populations with rho=1 and (ux,uy)=(0,0)
feq=zeros(NPOP);
f1=zeros(NPOP,NX,NY);
f2=zeros(NPOP,NX,NY);
forcepop=zeros(NPOP);
for y=1:NY
    for x=1:NX
        dense=rho(x,y);
        vx=ux(x,y);
        vy=uy(x,y);
        for k=1:NPOP
            feq(k)=weights(k)*(dense+(3*(vx*cx(k)+vy*cy(k)) ...
                +9/2*(cx(k)*vx+cy(k)*vy)^2-3/2*(vx^2+vy^2)));
            f1(k,x,y)=feq(k);
            f2(k,x,y)=feq(k);
        end
    end
end

% Main algorithm
for counter=1:NSTEPS

    % Macroscopic parameters computed through velocity moments of
    % populations f1
     for y=1:NY
        for x=1:NX
            
            dense=0;
            vx=0;
            vy=0;
            for k=1:NPOP
                dense=dense+f1(k,x,y);
                vx=vx+cx(k)*f1(k,x,y);
                vy=vy+cy(k)*f1(k,x,y);
            end
            
            rho(x,y)=dense;
            ux(x,y)=vx;
            uy(x,y)=vy;


            for k=1:NPOP

                % Compute the populations equilibrium value
                feq(k)=weights(k).*(dense+3*(vx*cx(k)+vy*cy(k)) ...
                    +9/2*((cx(k)*cx(k)-1/3)*vx*vx+2*cx(k)*cy(k)*vx*vy+(cy(k)*cy(k)-1/3)*vy*vy));
                
                % Compute external forcing term
                forcepop(k)=weights(k).*3.*(cx(k).*forcex+cy(k).*forcey);


                % Collision step
                f1(k,x,y)=f1(k,x,y)*(1-omega)+feq(k)*omega+forcepop(k);
                
                % Streaming step               
                newx=1+mod(x-1+cx(k)+NX,NX);
                newy=1+mod(y-1+cy(k)+NY,NY);
                f2(k,newx,newy)=f1(k,x,y);
            end
        end
     end


    % Bounceback Boundary Conditions
    
    for y=1:NY
        for x=1:NX
            if y==1 % Bottom wall

             % Introduce here the BB boundary conditions!   
                
            end
            if y==NY % Top wall

             % Introduce here the BB boundary conditions!   
             
            end
        end
    end

    % Assign new state f1, i.e. f(t+1) to previous state f2, i.e. f(t)
    f1=f2;
end

ux_plot=zeros(NX,NY+2);
ux_plot(:,2:NY+1)=ux;
% Analytical solution
y_plot=[y_bottom,1:NY,y_top];
ux_analy=-1/(2*kvisc).*forcex.*(y_plot-y_bottom).*(y_plot-y_top);


% % Compare LBM and Analytical solution velocity profiles
% 
% figure('color',[1 1 1])
% hold on
% plot(y_plot./(y_top),ux_analy./umax,'ko--');
% xlabel('y/y_{top}');
% ylabel('u/umax');
% plot(y_plot./(y_top),ux_plot(round(NX/2),:)./umax,'rs-.');
% legend('ux analy','ux LBM');
% axis  tight
% box on


% Calculation of L2 error
sum_num=0;
sum_denom=0;
for y=1:NY+2
    for x=1:NX
        sum_num=sum_num+(ux_plot(x,y)-ux_analy(y)).^2;
        sum_denom=sum_denom+ux_analy(y).^2;
    end
end

error=sqrt((sum_num)/(sum_denom));

disp(['L2 relative error = ',num2str(error)]);

toc % Stop time counter

