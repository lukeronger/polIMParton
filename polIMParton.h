#ifndef POLIMPARTON_H_
#define POLIMPARTON_H_ 1

class polIMParton
//this class is for users
{
public:
	polIMParton(unsigned int Z_temp=1, unsigned int A_temp=1);               //constructor function, the first parameter is Z, the other is A for a nuclei
	virtual void setDataSet(int);                            //Choose a data set, 1 is for set A and 2 is for set B
	virtual double getPolPDF(int, double x, double Q2) const;   //user function to get parton distribution functions, see ReadMe.txt for details       
        virtual double getXUV(double x, double Q2) const;        //return x(u -ubaar)
        virtual double getXDV(double x, double Q2) const;        //return x(d - dbar)
        virtual double getXUSea(double x, double Q2) const;      //return 2x*ubar
        virtual double getXDSea(double x, double Q2) const;      //return 2x*dbar
        virtual double getXSSea(double x, double Q2) const;      //return 2x*sbar
        virtual double getXCSea(double x, double Q2) const;      //return 2x*cbar
        virtual double getXGluon(double x, double Q2) const;     //return x*gluon
	virtual double getPolPDFError(int, double x, double Q2) const;   //user function to get the error of parton distribution function, see ReadMe.txt for details       
        virtual double getXUVError(double x, double Q2) const;        //return error of x(u -ubaar)
        virtual double getXDVError(double x, double Q2) const;        //return error of x(d - dbar)
        virtual double getXUSeaError(double x, double Q2) const;      //return error of 2x*ubar
        virtual double getXDSeaError(double x, double Q2) const;      //return error of 2x*dbar
        virtual double getXSSeaError(double x, double Q2) const;      //return error of 2x*sbar
        virtual double getXCSeaError(double x, double Q2) const;      //return error of 2x*cbar
        virtual double getXGluonError(double x, double Q2) const;     //return error of x*gluon
	virtual ~polIMParton(void);                                 //deconstructor function

private:
	unsigned int xMax;        //grid points number for variable x
	unsigned int Q2Max;       //grid points number for variable Q^2
	unsigned int flavorMax;   //flavor number in the grid data
	double lnxstep;           //step length of ln(x) for the grid data
	double lnQ2step;          //step length of ln(Q^2) for the grid data
	unsigned int Z;           //atomic number for a nuclei
	unsigned int A;           //mass number for a nuclei

        double * grid;            //grid data array
        double * griderr;         //grid error data array
        double * gridA;           //data array storing data set A
        double * gridAError;      //data array storing error data set A

	double fitLinear(double x, double* px, double* pf) const;       //linear interpolation function
	double fitQuadratic(double x, double* px, double* pf) const;    //quadratic interpolation function
        double getPDFType(int, double x, double Q2) const;
        double getPDFErrorType(int, double x, double Q2) const;
  
};

#endif
