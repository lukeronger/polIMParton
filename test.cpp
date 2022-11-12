#include<iostream>
#include<fstream>
extern "C"{
#include<math.h>
#include<time.h>
#include<unistd.h>
#include<sys/time.h>
}
#include "polIMParton.h"
using namespace std;


polIMParton proton(1,1);              // Z=1 and A=1, which is proton.
polIMParton neutron(0,1);             // Z=0 and A=1, which is neutron.




int main()
{
    //using data set A
    proton.setDataSet(1);             //1 for data set A and 2 data set B
    neutron.setDataSet(1);             //1 for data set A and 2 data set B

    struct timeval tstart,tend;
    gettimeofday(&tstart,NULL);       //get the starting time



    double lnxstep=(log(1)-log(0.1))/100;       //step length of log(x) for output data
    //output file for polarized xuv of proton at Q^2 = 1 GeV^2
    ofstream outdata1("proton_Q2_1GeV2_xDeltauv.txt");
    //output file for polarized xgluon of neutron at Q^2 = 100 GeV^2
    ofstream outdata2("neutron_Q2_100GeV2_xDeltagluon.txt");
    //get the polarized parton distribution functions and errors
    for(int i=0;i<700;i++)
    {
        double x = exp(log(1e-7)+i*lnxstep);
        outdata1<<x<<" "<<  x*(proton.getPolPDF(1, x,1.0)-proton.getPolPDF(-1, x,1.0));
	outdata1<<" "<<0<<" "<< x*(proton.getPolPDFError(2, x,1.0)+proton.getPolPDFError(-2, x,1.0)) <<endl;
        outdata2<<x<<" "<<  x*(neutron.getPolPDF(0, x,100.0))  <<" "<<0<<" "<<x*(neutron.getPolPDFError(0, x,100.0))<<endl;
    }





    gettimeofday(&tend,NULL);         //get the ending time
    double timeuse=1000000*(tend.tv_sec-tstart.tv_sec)+(tend.tv_usec-tstart.tv_usec);
    //display the total program running time
    cout<<"    Runtime : "<<(timeuse/1000.0)<<" ms."<<endl;


}



