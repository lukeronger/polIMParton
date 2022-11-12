#include<iostream>
#include<fstream>
#include<string>
extern "C"{
#include<math.h>
}
#include "polIMParton.h"
using namespace std;


//a method used to choose a data set
void polIMParton::setDataSet(int dataset)
{
        if(dataset==1)
        {
                grid = gridA;
		griderr = gridAError;
                cout<<"    Using data set A."<<endl;
        }
        else
        {
                cout<<"!!->Unknown data set."<<endl;
                cout<<"!!->Data set should be 1 only, so far"<<endl;
        }
}

//return the parton distributions of different kinds at x and Q^2
double polIMParton::getPolPDF(int Iparton, double x, double Q2) const
{
        if(Iparton==-4 || Iparton==4)return getXCSea(x,Q2)/2.0/x;
        else if(Iparton==-3 || Iparton==3)return getXSSea(x,Q2)/2.0/x;
        else if(Iparton==-2)return getXDSea(x,Q2)/2.0/x;
        else if(Iparton==2)return getXDSea(x,Q2)/2.0/x+getXDV(x,Q2)/x;
        else if(Iparton==-1)return getXUSea(x,Q2)/2.0/x;
        else if(Iparton==1)return getXUSea(x,Q2)/2.0/x+getXUV(x,Q2)/x;
        else if(Iparton==0)return getXGluon(x,Q2)/x;
        else
        {
                cout<<"!!->Unknown Iparton type."<<" (Iparton = "<<Iparton<<"?)"<<endl;
                cout<<"!!->Iparton should be one of these: [-4,-3, ... 3, 4]"<<endl;
                return 0;
        }
}
//return the parton distribution errors of different kinds at x and Q^2
double polIMParton::getPolPDFError(int Iparton, double x, double Q2) const
{
        if(Iparton==-4 || Iparton==4)return getXCSeaError(x,Q2)/2.0/x;
        else if(Iparton==-3 || Iparton==3)return getXSSeaError(x,Q2)/2.0/x;
        else if(Iparton==-2)return getXDSeaError(x,Q2)/2.0/x;
        else if(Iparton==2)return getXDSeaError(x,Q2)/2.0/x+getXDVError(x,Q2)/x;
        else if(Iparton==-1)return getXUSeaError(x,Q2)/2.0/x;
        else if(Iparton==1)return getXUSeaError(x,Q2)/2.0/x+getXUVError(x,Q2)/x;
        else if(Iparton==0)return getXGluonError(x,Q2)/x;
        else
        {
                cout<<"!!->Unknown Iparton type."<<" (Iparton = "<<Iparton<<"?)"<<endl;
                cout<<"!!->Iparton should be one of these: [-4,-3, ... 3, 4]"<<endl;
                return 0;
        }
}

//the constructor and initialization
polIMParton::polIMParton(unsigned int Z_temp, unsigned int A_temp)
:Z(Z_temp),A(A_temp)            //Z and A are parameters for a nuclei
{
	cout<<"    polIMParton version - 1.0"<<endl;
	cout<<"    Polarized Parton Distributions of Nucleons"<<endl;
	char filename[50];
	ifstream datain;
	double x, Q2;
	unsigned int i, j;
	xMax=61;
	Q2Max=32;
	flavorMax=7;
	lnxstep=log(10)/((xMax-1)/6);
	lnQ2step=log(2.0);
	gridA = new double[Q2Max*xMax*flavorMax];
	gridAError = new double[Q2Max*xMax*flavorMax];
	//read grid data for interpolation
	//reading data set A
	sprintf(filename,"grid_%d_%d_pol_SetA.dat",1,1);
	cout<<"    Loading "<<filename<<endl;
	datain.open(filename);
	if(!datain.good())cout<<"!!->Error while opening "<<filename<<"!\n!!->grid data file not exist?"<<endl;
	else
	for(i=0;i<Q2Max;i++)
	{
		for(j=0;j<xMax;j++)
		{
			datain>>Q2>>x>>(*(gridA+(xMax*i+j)*7))>>(*(gridA+(xMax*i+j)*7+1))>>(*(gridA+(xMax*i+j)*7+2))>>(*(gridA+(xMax*i+j)*7+3))>>(*(gridA+(xMax*i+j)*7+4))>>(*(gridA+(xMax*i+j)*7+5))>>(*(gridA+(xMax*i+j)*7+6));
		}
	}
	datain.close();
	//reading error data set A
	sprintf(filename,"grid_%d_%d_pol_err_SetA.dat",1,1);
	cout<<"    Loading "<<filename<<endl;
	datain.open(filename);
	if(!datain.good())cout<<"!!->Error while opening "<<filename<<"!\n!!->grid error data file not exist?"<<endl;
	else
	for(i=0;i<Q2Max;i++)
	{
		for(j=0;j<xMax;j++)
		{
			datain>>Q2>>x>>(*(gridAError+(xMax*i+j)*7))>>(*(gridAError+(xMax*i+j)*7+1))>>(*(gridAError+(xMax*i+j)*7+2))>>(*(gridAError+(xMax*i+j)*7+3))>>(*(gridAError+(xMax*i+j)*7+4))>>(*(gridAError+(xMax*i+j)*7+5))>>(*(gridAError+(xMax*i+j)*7+6));
		}
	}
	datain.close();
	//the default is set A
	grid = gridA;
	griderr = gridAError;
}

//the deconstructor
polIMParton::~polIMParton(void)
{
	delete[] gridA;
        delete[] gridAError;

}

//a method which returns xuv
double polIMParton::getXUV(double x, double Q2) const
{
	if(Z==0)return getPDFType(2,x,Q2);
	else return getPDFType(1,x,Q2); 
}
//a method which returns xdv
double polIMParton::getXDV(double x, double Q2) const
{
	if(Z==0)return getPDFType(1,x,Q2);
	else return getPDFType(2,x,Q2);
}
//a method which returns xusea
double polIMParton::getXUSea(double x, double Q2) const
{
	if(Z==0)return getPDFType(4,x,Q2);
	else return getPDFType(3,x,Q2);
}
//a method which returns xdsea
double polIMParton::getXDSea(double x, double Q2) const
{
	if(Z==0)return getPDFType(3,x,Q2);
	else return getPDFType(4,x,Q2);
}
//a method which returns xssea
double polIMParton::getXSSea(double x, double Q2) const
{
        return getPDFType(5,x,Q2);
}
//a method which returns xcsea
double polIMParton::getXCSea(double x, double Q2) const
{
        return getPDFType(6,x,Q2);
}
//a method which returns xgluon
double polIMParton::getXGluon(double x, double Q2) const
{
        return getPDFType(0,x,Q2);
}



//a method which returns the error of xuv
double polIMParton::getXUVError(double x, double Q2) const
{
	if(Z==0)return getPDFErrorType(2,x,Q2);
	else return getPDFErrorType(1,x,Q2); 
}
//a method which returns the error of xdv
double polIMParton::getXDVError(double x, double Q2) const
{
	if(Z==0)return getPDFErrorType(1,x,Q2);
	else return getPDFErrorType(2,x,Q2);
}
//a method which returns the error of xusea
double polIMParton::getXUSeaError(double x, double Q2) const
{
	if(Z==0)return getPDFErrorType(4,x,Q2);
	else return getPDFErrorType(3,x,Q2);
}
//a method which returns the error of xdsea
double polIMParton::getXDSeaError(double x, double Q2) const
{
	if(Z==0)return getPDFErrorType(3,x,Q2);
	else return getPDFErrorType(4,x,Q2);
}
//a method which returns the error of xssea
double polIMParton::getXSSeaError(double x, double Q2) const
{
        return getPDFErrorType(5,x,Q2);
}
//a method which returns the error of xcsea
double polIMParton::getXCSeaError(double x, double Q2) const
{
        return getPDFErrorType(6,x,Q2);
}
//a method which returns the error of xgluon
double polIMParton::getXGluonError(double x, double Q2) const
{
        return getPDFErrorType(0,x,Q2);
}




//a method which returns different types of distributions
double polIMParton::getPDFType(int Iparton, double x, double Q2) const
{
	if(Iparton<0 || Iparton>6)
	{
		cout<<"!!->Wrong Iparton input for getPDFType(int Iparton, double x, double Q2)."<<endl;
		return 0;
	}
	else
	{
		double lnx, lnQ2;
		int i=(int)(lnx=log(x*1e6)/lnxstep);
		int j=(int)(lnQ2=log(Q2*8)/lnQ2step);
		double g0[3], g1[3], g2[3], g[3]={0};
		if(i<0)i=0;
		if(i>(int)(xMax-3))i=xMax-3;
		if(j<0)j=0;
		if(j>29)j=29;
		//avoid log(1-x) calculation in below algorithm
		if(x>0.9999)return 0.0;
		//if x<2e-6, do linear interpolation or extrapolation
		else if(x<2e-6 || x>0.55)
		{
			double vx[2]={exp(log(1e-6)+i*lnxstep),  exp(log(1e-6)+i*lnxstep+lnxstep)};
			g0[0]=grid[(xMax*j+i)*7+Iparton];
			g0[1]=grid[(xMax*j+i+1)*7+Iparton];
			j++;
			g1[0]=grid[(xMax*j+i)*7+Iparton];
			g1[1]=grid[(xMax*j+i+1)*7+Iparton];
			j++;
			g2[0]=grid[(xMax*j+i)*7+Iparton];
			g2[1]=grid[(xMax*j+i+1)*7+Iparton];
			g[0]=fitLinear(x,vx,g0);
			g[1]=fitLinear(x,vx,g1);
			g[2]=fitLinear(x,vx,g2);
		}
		//we use quadratic interpolation method for other situations
		else
		{
			double vlnx[3]={(double)i,(double)(i+1),(double)(i+2)};
			g0[0]=grid[(xMax*j+i)*7+Iparton];
			g0[1]=grid[(xMax*j+i+1)*7+Iparton];
			g0[2]=grid[(xMax*j+i+2)*7+Iparton];
			j++;
			g1[0]=grid[(xMax*j+i)*7+Iparton];
			g1[1]=grid[(xMax*j+i+1)*7+Iparton];
			g1[2]=grid[(xMax*j+i+2)*7+Iparton];
			j++;
			g2[0]=grid[(xMax*j+i)*7+Iparton];
			g2[1]=grid[(xMax*j+i+1)*7+Iparton];
			g2[2]=grid[(xMax*j+i+2)*7+Iparton];
			g[0]=fitQuadratic(lnx,vlnx,g0);
			g[1]=fitQuadratic(lnx,vlnx,g1);
			g[2]=fitQuadratic(lnx,vlnx,g2);
		}
		//if Q2>1, we do the interpolation to the variable ln(Q^2)
		if(Q2>1)
		{
			double vlnQ2[3]={(double)(j-2),(double)(j-1),(double)j};
			return fitQuadratic(lnQ2,vlnQ2,g);
		}
		//if Q2<1, we do the interpolation to the variable Q^2
		else 
		{
			double vQ2[3]={0.125*pow(2,j-2),0.125*pow(2,j-1),0.125*pow(2,j)};
			return fitQuadratic(Q2,vQ2,g);
		}	
	}
}
//a method which returns different types of distribution errors
double polIMParton::getPDFErrorType(int Iparton, double x, double Q2) const
{
	if(Iparton<0 || Iparton>6)
	{
		cout<<"!!->Wrong Iparton input for getPDFErrorType(int Iparton, double x, double Q2)."<<endl;
		return 0;
	}
	else
	{
		double lnx, lnQ2;
		int i=(int)(lnx=log(x*1e6)/lnxstep);
		int j=(int)(lnQ2=log(Q2*8)/lnQ2step);
		double g0[3], g1[3], g2[3], g[3]={0};
		if(i<0)i=0;
		if(i>(int)(xMax-3))i=xMax-3;
		if(j<0)j=0;
		if(j>29)j=29;
		//avoid log(1-x) calculation in below algorithm
		if(x>0.9999)return 0.0;
		//if x<2e-6, do linear interpolation or extrapolation
		else if(x<2e-6 || x>0.55)
		{
			double vx[2]={exp(log(1e-6)+i*lnxstep),  exp(log(1e-6)+i*lnxstep+lnxstep)};
			g0[0]=griderr[(xMax*j+i)*7+Iparton];
			g0[1]=griderr[(xMax*j+i+1)*7+Iparton];
			j++;
			g1[0]=griderr[(xMax*j+i)*7+Iparton];
			g1[1]=griderr[(xMax*j+i+1)*7+Iparton];
			j++;
			g2[0]=griderr[(xMax*j+i)*7+Iparton];
			g2[1]=griderr[(xMax*j+i+1)*7+Iparton];
			g[0]=fitLinear(x,vx,g0);
			g[1]=fitLinear(x,vx,g1);
			g[2]=fitLinear(x,vx,g2);
		}
		//we use quadratic interpolation method for other situations
		else
		{
			double vlnx[3]={(double)i,(double)(i+1),(double)(i+2)};
			g0[0]=griderr[(xMax*j+i)*7+Iparton];
			g0[1]=griderr[(xMax*j+i+1)*7+Iparton];
			g0[2]=griderr[(xMax*j+i+2)*7+Iparton];
			j++;
			g1[0]=griderr[(xMax*j+i)*7+Iparton];
			g1[1]=griderr[(xMax*j+i+1)*7+Iparton];
			g1[2]=griderr[(xMax*j+i+2)*7+Iparton];
			j++;
			g2[0]=griderr[(xMax*j+i)*7+Iparton];
			g2[1]=griderr[(xMax*j+i+1)*7+Iparton];
			g2[2]=griderr[(xMax*j+i+2)*7+Iparton];
			g[0]=fitQuadratic(lnx,vlnx,g0);
			g[1]=fitQuadratic(lnx,vlnx,g1);
			g[2]=fitQuadratic(lnx,vlnx,g2);
		}
		//if Q2>1, we do the interpolation to the variable ln(Q^2)
		if(Q2>1)
		{
			double vlnQ2[3]={(double)(j-2),(double)(j-1),(double)j};
			return fitQuadratic(lnQ2,vlnQ2,g);
		}
		//if Q2<1, we do the interpolation to the variable Q^2
		else 
		{
			double vQ2[3]={0.125*pow(2,j-2),0.125*pow(2,j-1),0.125*pow(2,j)};
			return fitQuadratic(Q2,vQ2,g);
		}	
	}
}


//linear interpolation method
double polIMParton::fitQuadratic(double x, double* px, double* pf) const
{
	double f01=(pf[1]-pf[0])/(px[1]-px[0]);
	double f12=(pf[2]-pf[1])/(px[2]-px[1]);
	double f012=(f12-f01)/(px[2]-px[0]);
	return pf[0]+f01*(x-px[0])+f012*(x-px[0])*(x-px[1]);
}

//quadratic interpolation method
double polIMParton::fitLinear(double x, double* px, double* pf) const
{
	double f01=(pf[1]-pf[0])/(px[1]-px[0]);
	return pf[0]+f01*(x-px[0]);
}




