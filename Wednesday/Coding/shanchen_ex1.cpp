#include <iostream>
#include <cmath>
#include <time.h>
#include <fstream>
#include <sstream>
#include <string>
#include <cstdlib>
#include <algorithm>

//Parameters of the grid
const int nx=128;
const int ny=128;
const int npop=9;
const int nsteps=7000;
const int noutput=100;

//Parameters of the LBM
const int cx[]={0,1,0,-1,0,1,-1,-1,1};
const int cy[]={0,0,1,0,-1,1,1,-1,-1};
const double weights[]={4.0/9.0,1.0/9.0,1.0/9.0,1.0/9.0,1.0/9.0,1.0/36.0,1.0/36.0,1.0/36.0,1.0/36.0};
double tau=1.0;

//Parameters of the Shan-Chen model
double rhol=1.95;
double rhog=0.15;
int radius=30;
double g=-5.0;

//Arrays
double rho[nx][ny];
double u1[nx][ny];
double u2[nx][ny];
double f[nx][ny][npop];
double f2[nx][ny][npop];
double feq[npop];

void writeFile(std::string name,double (*fun)[nx][ny])
{
	std::ofstream fout(name.c_str());
	for(int counterX=0;counterX<nx;counterX++)
	{	
		for(int counterY=0;counterY<ny;counterY++)
			fout<<(*fun)[counterX][counterY]<<" ";
		fout<<"\n";
	}		
	fout<<std::endl;
}

void calculateDensities()
{
	double rho_gas=0.0;
	double rho_liq=0.0;
	
	for(int counterX=0;counterX<10;counterX++)
		for(int counterY=0;counterY<10;counterY++)
			rho_gas+=rho[counterX][counterY];
	
	for(int counterX=nx/2-5;counterX<nx/2+5;counterX++)
		for(int counterY=ny/2-5;counterY<ny/2+5;counterY++)
			rho_liq+=rho[counterX][counterY];

	rho_gas=rho_gas/100.0;
	rho_liq=rho_liq/100.0;
	
	double pressure_liq,pressure_gas;
	pressure_liq=rho_liq/3.0+g/6.0*(1.0-exp(-rho_liq))*(1.0-exp(-rho_liq));
	pressure_gas=rho_gas/3.0+g/6.0*(1.0-exp(-rho_gas))*(1.0-exp(-rho_gas));
	std::cout<<"Calculated densities are: Rho_gas="<<rho_gas<<" Rho_liq="<<rho_liq<<"\n";
	std::cout<<"Delta Pressure is:"<<pressure_liq-pressure_gas<<"\n";
}


int main(int argc, char** argv)
{
	
	if (argc==4)
	{
		g=atof(argv[1]);
		tau=atof(argv[2]);
		radius=atoi(argv[3]);
	}
	else if (argc==6)
	{
		g=atof(argv[1]);
		tau=atof(argv[2]);
		radius=atoi(argv[3]);
		rhog=atof(argv[4]);
		rhol=atof(argv[5]);
	}

	
	//Initialization
	for(int counterX=0;counterX<nx;counterX++)
		for(int counterY=0;counterY<ny;counterY++)
		{	
			if ((counterY-ny/2.0)*(counterY-ny/2.0)+(counterX-nx/2.0)*(counterX-nx/2.0)<=radius*radius)
			{
				rho[counterX][counterY]=rhol;
			}
			else 
				rho[counterX][counterY]=rhog;

			double dense,v1,v2;
		
			dense=rho[counterX][counterY];
			v1=v2=u1[counterX][counterX]=u2[counterX][counterY]=0.0;
			double usq = v1*v1 + v2*v2;
			feq[0] = 4.0/9.0 * dense * (1.0 - 1.5 * usq); 
			feq[1] = 1.0/9.0 * dense * (1.0 + 3*v1 + 4.5*v1*v1 - 1.5*usq); 
			feq[2] = 1.0/9.0 * dense * (1.0 + 3*v2 + 4.5*v2*v2 - 1.5*usq); 
			feq[3] = 1.0/9.0 * dense * (1.0 - 3*v1 + 4.5*v1*v1 - 1.5*usq); 
			feq[4] = 1.0/9.0 * dense * (1.0 - 3*v2 + 4.5*v2*v2 - 1.5*usq); 
			feq[5] = 1.0/36.0 * dense * (1.0 + 3*(v1 + v2) + 4.5*(v1 + v2)*(v1 + v2) - 1.5*usq); 
			feq[6] = 1.0/36.0 * dense * (1.0 + 3*(-v1 + v2) + 4.5*(-v1 + v2)*(-v1 + v2) - 1.5*usq);
			feq[7] = 1.0/36.0 * dense * (1.0 + 3*(-v1 - v2) + 4.5*(v1 + v2)*(v1 + v2) - 1.5*usq); 
			feq[8] = 1.0/36.0 * dense * (1.0 + 3*(v1 - v2) + 4.5*(v1 - v2)*(v1 -v2) - 1.5*usq); 
			for (int k=0; k<npop; k++) {
				f[counterX][counterY][k]=feq[k];
				f2[counterX][counterY][k]=feq[k];
			}
		}
	
	time_t start, finish;
	start = time(NULL);
	
	//Main loop
	for (int timecounter=0; timecounter<=nsteps;timecounter++) 
	{
		
		//Calculation of the density field
		for(int counterX=0;counterX<nx;counterX++)
			for(int counterY=0;counterY<ny;counterY++)
			{
				rho[counterX][counterY]=0; 
				for (int k=0; k<9; k++ )
				{			
					rho[counterX][counterY]+=f[counterX][counterY][k]; 
				}		
				
			}

		double rho_temp[nx][ny];
		
		//Collision and streaming
		for(int counterX=0;counterX<nx;counterX++)
			for(int counterY=0;counterY<ny;counterY++)
			{

			
				double dense,v1,v2;
			
				dense=rho[counterX][counterY];

				double fx=0.0;
				double fy=0.0;

				for(int k=0;k<9;k++)
				{
					int iX2=(counterX+cx[k]+nx) % nx; 
					int iY2=(counterY+cy[k]+ny) % ny;
					fx+=weights[k]*cx[k]*(1.0-exp(-rho[iX2][iY2]));
					fy+=weights[k]*cy[k]*(1.0-exp(-rho[iX2][iY2]));
				}
			
				fx=-g*(1.0-exp(-dense))*fx;
				fy=-g*(1.0-exp(-dense))*fy;
			
				//v1=u1[i]=(f[9*i+1]-f[9*i+3]+f[9*i+5]-f[9*i+6]-f[9*i+7]+f[9*i+8])/dense+fx/(tau*dense); 
				//v2=u2[i]=(f[9*i+2]-f[9*i+4]+f[9*i+5]+f[9*i+6]-f[9*i+7]-f[9*i+8])/dense+fy/(tau*dense); 

				v1=u1[counterX][counterY]=(f[counterX][counterY][1]-f[counterX][counterY][3]
					+f[counterX][counterY][5]-f[counterX][counterY][6]-f[counterX][counterY][7]
					+f[counterX][counterY][8])/dense+fx/(2.0*dense); 
				v2=u2[counterX][counterY]=(f[counterX][counterY][2]-f[counterX][counterY][4]
					+f[counterX][counterY][5]+f[counterX][counterY][6]-f[counterX][counterY][7]
					-f[counterX][counterY][8])/dense+fy/(2.0*dense); 
			
				double fpop[9];
				for(int k=0;k<9;k++)
					fpop[k]=weights[k]*(1.0-0.5/tau)*((3.0*(cx[k]-v1)+9.0*cx[k]*(cx[k]*v1+cy[k]*v2))*fx
		               +(3.0*(cy[k]-v2)+9.0*cy[k]*(cx[k]*v1+cy[k]*v2))*fy);
			
			
				float usq = v1*v1 + v2*v2;	
			
				feq[0] = 4.0/9.0 * dense * (1.0 - 1.5 * usq); 
				feq[1] = 1.0/9.0 * dense * (1.0 + 3*v1 + 4.5*v1*v1 - 1.5*usq); 
				feq[2] = 1.0/9.0 * dense * (1.0 + 3*v2 + 4.5*v2*v2 - 1.5*usq); 
				feq[3] = 1.0/9.0 * dense * (1.0 - 3*v1 + 4.5*v1*v1 - 1.5*usq); 
				feq[4] = 1.0/9.0 * dense * (1.0 - 3*v2 + 4.5*v2*v2 - 1.5*usq); 
				feq[5] = 1.0/36.0 * dense * (1.0 + 3*(v1 + v2) + 4.5*(v1 + v2)*(v1 + v2) - 1.5*usq); 
				feq[6] = 1.0/36.0 * dense * (1.0 + 3*(-v1 + v2) + 4.5*(-v1 + v2)*(-v1 + v2) - 1.5*usq);
				feq[7] = 1.0/36.0 * dense * (1.0 + 3*(-v1 - v2) + 4.5*(v1 + v2)*(v1 + v2) - 1.5*usq);
				feq[8] = 1.0/36.0 * dense * (1.0 + 3*(v1 - v2) + 4.5*(v1 - v2)*(v1 -v2) - 1.5*usq);
			
				for(int k=0; k<9; k++) 
				{  
					f2[counterX][counterY][k]=f[counterX][counterY][k]-1.0/tau*(f[counterX][counterY][k]-feq[k])+fpop[k]; 
				}
				rho_temp[counterX][counterY]=dense;  

			}
		//if (timecounter==0)
		// 	writeFile("aftercollision.dat",&rho_temp);
		
		
		for(int counterX=0;counterX<nx;counterX++)
			for(int counterY=0;counterY<ny;counterY++)
			{
				for(int k=0;k<9;k++)
				{
					int iX2=(counterX+cx[k]+nx) % nx; 
					int iY2=(counterY+cy[k]+nx) % ny;
					f[iX2][iY2][k]=f2[counterX][counterY][k]; 
				}
	
			
			}
		
		
		//Preparation of the output
		if (timecounter%noutput==0)
		{
			std::stringstream imagestream;
			std::stringstream len;
			len<<timecounter;
		
			imagestream << "height"<<std::string(5-len.str().size(),'0')<<timecounter<<".dat";
			writeFile(imagestream.str(),&rho);	
			std::cout<<"Iteration is "<<timecounter<<"\n";
			calculateDensities();
		}
	}
	
	finish = time(NULL);
	
	std::cout<<"Overall time is "<<finish-start<<" sec"<<"\n";
	
    return 0;
}
