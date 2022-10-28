#include <iostream>
#include <fstream>
#include "twodproblem.h"

int n_interior_min;

// this  code is for testing smooth and non smooth test case for
// burger's equation in two dimensional case.
using namespace std;

void  compute_rate(const vector<double>& l1error, const double h )
{
	unsigned int n=l1error.size();

	std::cout<<"L1error"<<" "<<"Order "<<std::endl;
	for (unsigned int i=0;i<n;i++)
	{
		if(i==0)
      std::cout<<l1error[i]<<"  "<<"--"<<std::endl;
		else
      std::cout<<l1error[i]<<"  "<< log(l1error[i-1]/l1error[i])/log(2)<<std::endl;

	}
}
int main ()
{
   cout << "Starting reservoir problem ..." << endl;

   unsigned int n_cycles=1;
   std::ofstream  fo;
  // std::string filename="error.dat";

   fo.open("error.dat");
   fo<< "h" << " "<<" L1error "<<"  "<<" Order"<<std::endl;

   std::vector<double> L1error(n_cycles);

   	  for ( unsigned int j=0; j< n_cycles;j++)
   	  {
   		TwoDProblem twodproblem;

   		twodproblem.run (j,L1error[j]);

   	  }
     double h=0.0;
   	  compute_rate(L1error,h);
   	  for (unsigned int i=0; i< n_cycles;i++)
   	  {
   		  if(i==0)
   	      fo << h <<" "<<L1error[i]<<"  " <<"--"<<std::endl;
   		  else
   		  fo << h << " " << L1error[i]<<"  " <<log(L1error[i-1]/L1error[i])/log(2)<<std::endl;
   	  }

    fo.close();
}
