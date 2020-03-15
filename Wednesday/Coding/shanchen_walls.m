%% periodic flow my version
clear all;
clc;

% Macroscopic density and velocities
NX=128;
NY=128;
NPOP=9;
NSTEPS=1000;
NOUTPUT=10;

% Parameters of the Shan-Chen model
rho_crit=log(2)
rho_liq=1.95
rho_gas=0.15
radius=20
G=-5.0
gravity=0.0001
rho_wall=0.1

% Initialization 
rho=ones(NX,NY);
ux=zeros(NX,NY);
uy=zeros(NX,NY);

feq=zeros(NPOP);
f1=zeros(NPOP,NX,NY);
f2=zeros(NPOP,NX,NY);
forcex=zeros(NX,NY);
forcey=zeros(NY,NY);

% Parameters of the lattice
weights=[4/9 1/9 1/9 1/9 1/9 1/36 1/36 1/36 1/36];
cx=[0 1 0 -1 0 1 -1 -1 1];
cy=[0 0 1 0 -1 1 1 -1 -1];
omega=1.0;

%% Initialization section
for y=1:NY
    for x=1:NX
        %if abs(x-NX/2)<radius && y>NY-2*radius && y<NY
        if (x-NX/2)^2+(y-NY+radius+4)^2<=radius*radius
            rho(x,y)=rho_liq;
        else
            rho(x,y)=rho_gas;
        end
        
        ux(x,y)=0;
        uy(x,y)=0;
        vx=ux(x,y);
        vy=uy(x,y);
        
        for k=1:NPOP
            feq(k)=weights(k)*rho(x,y)*(1+3*(vx*cx(k)+vy*cy(k)) ...
               + 9/2*((cx(k)*cx(k)-1/3)*vx*vx+2*cx(k)*cy(k)*vx*vy+(cy(k)*cy(k)-1/3)*vy*vy));
            f1(k,x,y)=feq(k);
            f2(k,x,y)=feq(k);
        end
        
    end
end

%%Main loop without walls

for counter=1:300
    
    % Calculation of the macroscopic quantities
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
            ux(x,y)=vx/dense;
            uy(x,y)=vy/dense;
        end
    end    
    
    % Calculation of the force
    for y=1:NY
        for x=1:NX
            
            force_sum_x=0.0;
            force_sum_y=0.0;
            
            for k=2:NPOP
                newx=1+mod(x-1+cx(k)+NX,NX);
                newy=1+mod(y-1+cy(k)+NY,NY);
                psi=1-exp(-rho(newx,newy));
                force_sum_x=force_sum_x-G*weights(k)*psi*cx(k);
                force_sum_y=force_sum_y-G*weights(k)*psi*cy(k);
            end
            forcex(x,y)=(1-exp(-rho(x,y)))*force_sum_x;
            forcey(x,y)=(1-exp(-rho(x,y)))*force_sum_y;
            
        end
    end
    
    % main loop
    for y=1:NY
        for x=1:NX
             
            dense=rho(x,y);
            vx=ux(x,y)+0.5*forcex(x,y)/dense;
            vy=uy(x,y)+0.5*forcey(x,y)/dense;
            for k=1:NPOP
                feq(k)=weights(k)*dense*(1+3*(vx*cx(k)+vy*cy(k)) ...
                   +9/2*((cx(k)*cx(k)-1/3)*vx*vx+2*cx(k)*cy(k)*vx*vy+(cy(k)*cy(k)-1/3)*vy*vy));
               

                forcepop=weights(k)*(1-0.5*omega)*((3*(cx(k)-vx)+9*cx(k)*(cx(k)*vx+cy(k)*vy))*forcex(x,y)...
                   +(3*(cy(k)-vy)+9*cy(k)*(cx(k)*vx+cy(k)*vy))*forcey(x,y));
               
                newx=1+mod(x-1+cx(k)+NX,NX);
                newy=1+mod(y-1+cy(k)+NY,NY);
         
                f1(k,x,y)=f1(k,x,y)*(1-omega)+feq(k)*omega+forcepop;
                f2(k,newx,newy)=f1(k,x,y);
            end
           
        end
    end
    
    counter
    f1=f2;
end

for x=1:NX
    rho(x,1)=rho_wall;
    rho(x,NY)=rho_wall;
end

%% Main loop
counter_frame=1;
for counter=1:NSTEPS
    
    % Calculation of the macroscopic quantities
    for y=2:NY-1
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
            ux(x,y)=vx/dense;
            uy(x,y)=vy/dense;
        end
    end
    
    % Putting walls density
    for x=1:NX
        rho(x,1)=rho_wall;
        rho(x,NY)=rho_wall;
        ux(x,1)=0;
        uy(x,1)=0;
        ux(x,NY)=0;
        uy(x,NY)=0;
    end
    
    
    % Calculation of the force
    for y=2:NY-1
        for x=1:NX
            
            force_sum_x=0.0;
            force_sum_y=0.0;
            
            for k=2:NPOP
                newx=1+mod(x-1+cx(k)+NX,NX);
                newy=1+mod(y-1+cy(k)+NY,NY);
                psi=1-exp(-rho(newx,newy));
                force_sum_x=force_sum_x-G*weights(k)*psi*cx(k);
                force_sum_y=force_sum_y-G*weights(k)*psi*cy(k);
            end
            forcex(x,y)=(1-exp(-rho(x,y)))*force_sum_x;
            forcey(x,y)=(1-exp(-rho(x,y)))*force_sum_y+gravity;
            
        end
    end
    
    % main loop
    for y=2:NY-1
        for x=1:NX
             
            dense=rho(x,y);
            vx=ux(x,y)+0.5*forcex(x,y)/dense;
            vy=uy(x,y)+0.5*forcey(x,y)/dense;
            for k=1:NPOP
                feq(k)=weights(k)*dense*(1+3*(vx*cx(k)+vy*cy(k)) ...
                   +9/2*((cx(k)*cx(k)-1/3)*vx*vx+2*cx(k)*cy(k)*vx*vy+(cy(k)*cy(k)-1/3)*vy*vy));
               

                forcepop=weights(k)*(1-0.5*omega)*((3*(cx(k)-vx)+9*cx(k)*(cx(k)*vx+cy(k)*vy))*forcex(x,y)...
                   +(3*(cy(k)-vy)+9*cy(k)*(cx(k)*vx+cy(k)*vy))*forcey(x,y));
               
                newx=1+mod(x-1+cx(k)+NX,NX);
                newy=1+mod(y-1+cy(k)+NY,NY);
         
                f1(k,x,y)=f1(k,x,y)*(1-omega)+feq(k)*omega+forcepop;
                f2(k,newx,newy)=f1(k,x,y);
            end
           
        end
    end
    
    % prescribing BB populations and propagating them
    for x=1:NX
        % top wall 
        set=[2,5,6];
        opp=[4,7,8];
        for k=1:size(set)
            f2(opp(k),x,NY-1)=f1(set(k),x,NY-1);
        end
        % bottom wall
        set=[4,7,8];
        opp=[2,5,6];
        for k=1:size(set)
            f2(opp(k),x,2)=f1(set(k),x,2);
        end
    end
    
    f1=f2;
    counter

    if mod(counter,NOUTPUT)==0
        imagesc(rho');
        F(counter_frame) = getframe;
        counter_frame=counter_frame+1;
    end
end

movie(F,10)

disp('Rho_liq=')
disp(mean(mean(rho(NX/2-5:NX/2+5,NY/2-5:NY/2+5))))
disp('Rho_gas=')
disp(mean(mean(rho(1:10,1:10))))
