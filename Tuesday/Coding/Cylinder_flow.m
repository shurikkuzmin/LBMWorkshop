clear all
clc
tic
% flow around circular cylinder 

% Lattice parameters
weights=[4/9 1/9 1/9 1/9 1/9 1/36 1/36 1/36 1/36];
cx=[0 1 0 -1 0 1 -1 -1 1];
cy=[0 0 1 0 -1 1 1 -1 -1];
opp = [1 4 5 2 3 8 9 6 7];
% Numerical parameters
NX=50;  % Number of grids points along x
NY=50;  % Number of grid points along y
NPOP=9; % Number of populations used in velocity space discretization
NSTEPS=5000;    % Number of time steps/iterations

% Simulation parameters
y_bottom=1.5;
y_top=NY-0.5;


obst_r = (y_top+y_bottom)/20+1;  % radius of the cylinder
obst_x = NX/2+obst_r;   % position of the cylinder; (exact
obst_y = (y_top+y_bottom)/2+1;   % y-symmetry is avoided)


Re=100;  % Reynolds number
% omega=32/(20+sqrt(208));   % Relaxation frequency
omega=1/0.6;
kvisc=1/3*(1/omega-0.5); % Kinematic viscosity
% kvisc=1/sqrt(48);
% omega=1/(3*kvisc+0.5);
umax=Re*kvisc/((y_top-y_bottom)) ;% umax=0.001; % Mach number (can be understood as a CFL number)



forcex=8.*umax*kvisc./((y_top-y_bottom).^2);

% Analytical solution
y_plot=[y_bottom,1:NY,y_top];
ux_analy=-1/(2*kvisc).*forcex.*(y_plot-y_bottom).*(y_plot-y_top);

% Macroscopic parameters
rho0=1;
rho=ones(NX,NY);
ux=zeros(NX,NY);
uy=zeros(NX,NY);


fluid=ones(NX,NY);
for x=1:NX
    for y=1:NY
        if (x-obst_x)^2 + (y-obst_y)^2 <= obst_r^2;
            fluid(x,y)=0;
            
        end
    end
end



normal_x=zeros(NX-1,NY);
normal_y=zeros(NX,NY-1);
% Outer normal
for x=1:NX
    for y=1:NY
        if x<NX
        normal_x(x,y)=fluid(x+1,y)-fluid(x,y);
        end
        if y<NY
        normal_y(x,y)=fluid(x,y+1)-fluid(x,y);
        end
    end
end


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
            feq(k)=weights(k)*(dense+rho0*(3*(vx*cx(k)+vy*cy(k)) ...
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


                    % Collision step
                    f1(k,x,y)=f1(k,x,y)*(1-omega)+feq(k)*omega;

                    % Streaming step
                    newx=1+mod(x-1+cx(k)+NX,NX);
                    newy=1+mod(y-1+cy(k)+NY,NY);
                    f2(k,newx,newy)=f1(k,x,y);
                end
        end
    end

    % Boundaries
    %==================================================================

    % Zou He at Inlet and Outlet
    %==================================================================
    x=1;    % Inlet
    for y=1:NY

        ux(x,y)=ux_analy(y+1);
        uy(x,y)=0;
        rho(x,y)=ux(x,y)+(f2(1,x,y)+f2(3,x,y)+f2(5,x,y)+2*(f2(4,x,y)+f2(7,x,y)+f2(8,x,y)));

        f2(2,x,y)=f2(4,x,y)+2/3*ux(x,y);
        f2(6,x,y)=f2(8,x,y)+1/6*ux(x,y)+0.5*(f2(5,x,y)-f2(3,x,y))+1/2*uy(x,y);
        f2(9,x,y)=f2(7,x,y)+1/6*ux(x,y)-0.5*(f2(5,x,y)-f2(3,x,y))-1/2*uy(x,y);

    end

    x=NX;    % Outlet
    for y=1:NY

        ux(x,y)=ux(x-1,y);
        uy(x,y)=0;
        rho(x,y)=-ux(x,y)+(f2(1,x,y)+f2(3,x,y)+f2(5,x,y)+2*(f2(2,x,y)+f2(6,x,y)+f2(9,x,y)));

        f2(4,x,y)=f2(2,x,y)-2/3*ux(x,y);
        f2(8,x,y)=f2(6,x,y)-1/6*ux(x,y)-0.5*(f2(5,x,y)-f2(3,x,y))-1/2*uy(x,y);
        f2(7,x,y)=f2(9,x,y)-1/6*ux(x,y)+0.5*(f2(5,x,y)-f2(3,x,y))+1/2*uy(x,y);

    end


    % Bounceback Boundary Conditions at bottom and top walls and side walls
    %==================================================================
    for k=1:NPOP

        for x=2:NX-1
            for y=2:NY-1


                if normal_x(x,y)<0

                    if cx(k)<0

                        f2(k,x,y)=f1(opp(k),x,y);
                    end

 
                elseif normal_x(x-1,y)>0

                    if cx(k)>0

                        f2(k,x,y)=f1(opp(k),x,y);

                    end


                    if normal_y(x,y)<0

                        if cy(k)<0

                            f2(k,x,y)=f1(opp(k),x,y);
                        end

                        % Right wall
                    elseif normal_y(x,y-1)>0

                        if cy(k)>0

                            f2(k,x,y)=f1(opp(k),x,y);

                        end
                    end

                end
            end
        end
    end
           
%     contourf(ux);
%     drawnow

    % Assign new state f1, i.e. f(t+1) to previous state f2, i.e. f(t)
    f1=f2;
    

 
end


[x_plot,y_plot]=meshgrid(1:NX,1:NY);
ux_plot=zeros(NX,NY);
uy_plot=zeros(NX,NY);


for x=1:NX
    for y=1:NY
        if fluid(x,y)==1 
            ux_plot(x,y)=ux(x,y);
            uy_plot(x,y)=uy(x,y);
        else
            ux_plot(x,y)=NaN;
            uy_plot(x,y)=NaN;
        end
    end
end
toc % Stop time counter


% Plot velocity field
figure('color',[1 1 1])
set(gcf,'Position',[50 120 910 550]); clf;
quiver(x_plot,y_plot,uy_plot,ux_plot,'r');
hold on;
set(gca,'Color','k');
%Draw cylinder
theta=linspace(0,2*pi,51);
plot(obst_y+obst_r*cos(theta),obst_x+obst_r*sin(theta),'Color',[0.9 0.5 0.1],'Linewidth',2);
axis equal tight
xlabel('\it{y}');ylabel('\it{x}');
% Streamlines
sx=linspace(1,NX,15).^1;
[sx sy]=meshgrid(sx, 1);
verts=stream2(x_plot,y_plot,uy_plot,ux_plot,sx,sy);
streamline(verts);

% Particles animation
nparts=200;
streamparticles(verts, nparts,'Animate',10,'Framerate',50,'MarkerFaceColor','g','ParticleAlignment','off');



