%% periodic flow my version
clear all;
clc;

% Macroscopic density and velocities
NX=64;
NY=64;
NPOP=9;
NSTEPS=400;

rho0=1;
umax=0.001;

rho=ones(NX,NY);
ux=zeros(NX,NY);
uy=zeros(NX,NY);

uxinit=zeros(NX,NY);
uyinit=zeros(NX,NY);

feq=zeros(NPOP);
f1=zeros(NPOP,NX,NY);
f2=zeros(NPOP,NX,NY);

weights=[4/9 1/9 1/9 1/9 1/9 1/36 1/36 1/36 1/36];
cx=[0 1 0 -1 0 1 -1 -1 1];
cy=[0 0 1 0 -1 1 1 -1 -1];
omega=0.1;

for y=1:NY
    for x=1:NX
        rho(x,y)=rho0+3*0.25*umax^2*...
            (cos(4*pi*(x-1)/NX)-cos(4*pi*(y-1)/NY));
        ux(x,y)=umax*sin(2*pi*(x-1)/NX)*sin(2*pi*(y-1)/NY);
        uy(x,y)=umax*cos(2*pi*(x-1)/NX)*cos(2*pi*(y-1)/NY);
        uxinit(x,y)=umax*sin(2*pi*(x-1)/NX)*sin(2*pi*(y-1)/NY);
        uyinit(x,y)=umax*cos(2*pi*(x-1)/NX)*cos(2*pi*(y-1)/NY);
       
        vx=ux(x,y); 
        vy=uy(x,y); 
        for k=1:NPOP
            feq(k)=weights(k)*(rho(x,y)+3*(vx*cx(k)+vy*cy(k)) ...
               + 9/2*((cx(k)*cx(k)-1/3)*vx*vx+2*cx(k)*cy(k)*vx*vy+(cy(k)*cy(k)-1/3)*vy*vy));
            f1(k,x,y)=feq(k);
            f2(k,x,y)=feq(k);
        end
    end
end

track=zeros(1,NSTEPS);

%this value we need to obtain - it should be constant over the simulation
%umax*sin(2*pi*(NX/4-1)/NX)*sin(2*pi*(NY/4-1)/NY)

for counter=1:NSTEPS
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
            vx=vx/dense;
            vy=vy/dense;
            
            % Collision
            %forcex=tempx(x,y);
            %forcey=tempy(x,y);
            forcex=(2*pi/NX)^2*2*1/3*(1/omega-1/2)*umax*sin(2*pi*(x-1)/NX)*sin(2*pi*(y-1)/NY);
            forcey=(2*pi/NX)^2*2*1/3*(1/omega-1/2)*umax*cos(2*pi*(x-1)/NX)*cos(2*pi*(y-1)/NY);
            
             
            rho(x,y)=dense;
            ux(x,y)=vx+forcex/(2.0*dense);
            uy(x,y)=vy+forcey/(2.0*dense);
 
            vx=vx+forcex/2.0;
            vy=vy+forcey/2.0;
 
            
            for k=1:NPOP
                feq(k)=weights(k)*dense*(1+3*(vx*cx(k)+vy*cy(k)) ...
                   +9/2*((cx(k)*cx(k)-1/3)*vx*vx+2*cx(k)*cy(k)*vx*vy+(cy(k)*cy(k)-1/3)*vy*vy));
                
                newx=1+mod(x-1+cx(k)+NX,NX);
                newy=1+mod(y-1+cy(k)+NY,NY);
         
                
                forcepop=weights(k)*(1-0.5*omega)*(forcex*(3.0*(cx(k)-vx)+9.0*cx(k)*(cx(k)*vx+cy(k)*vy))+...
                    forcey*(3.0*(cy(k)-vy)+9.0*cy(k)*(cx(k)*vx+cy(k)*vy)));
                f1(k,x,y)=f1(k,x,y)*(1-omega)+feq(k)*omega+forcepop;
                f2(k,newx,newy)=f1(k,x,y);
            end
           
        end
    end
    track(counter)=ux(NX/4+1,NY/4+1);
    %decay(counter)=umax*exp(-1/3*(1/omega-0.5)*counter*4*pi/NX);
    f1=f2;
end
figure(1)
surf(ux-uxinit)
figure(2)
plot(1:NSTEPS,track)

max(max(ux))
