#include"rotate.h"
#include"matrix_functions.h"
#include"printing_functions.h"
#include<omp.h>
#include"qmcinfo.h"
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multimin.h>

double cost(std::vector<RMatrix> &onerdm, double alpha1, double alpha2, 
	                                  double alpha3, double alpha4, bool print_on)
{
   std::vector<RMatrix> in_one_rdm = onerdm;
   int nstates=in_one_rdm.size();
   std::vector<RMatrix> out_one_rdm(nstates);
   RMatrix unitary;
   int nsites=12;
   unitary.resize(nsites,nsites);
   //cout<<"Setting unitary in cost"<<endl;
   for (int i=0;i<nsites;i++)
   {
	for (int j=0;j<nsites;j++)
	{
   		unitary(i,j)=0.0; 
	}
   }
   for (int i=0;i<4;i++) 
   {unitary(i,i)=sqrt(abs(1.0-(4.0*alpha1*alpha1)-(2.0*alpha2*alpha2)-(alpha4*alpha4)-(4.0*alpha3*alpha3)));}

//===================================================================================
/* 0     1       2      3      4      5      6      7     8      9      10      11
+0.913 +0.081 +0.006 +0.003 +0.200 +0.199 -0.201 -0.200 -0.028 +0.011 +0.012 -0.032 
+0.006 +0.913 -0.036 +0.088 -0.015 -0.200 +0.197 +0.035 -0.013 +0.015 +0.202 -0.198 
+0.083 +0.018 +0.913 +0.008 -0.200 +0.016 +0.001 +0.199 +0.199 -0.200 +0.055 -0.012 
-0.006 +0.006 +0.086 +0.913 -0.012 +0.041 -0.006 +0.021 -0.200 +0.198 -0.200 +0.202 
-0.199 +0.032 +0.198 -0.006 +0.949 +0.071 +0.046 +0.069 +0.070 +0.025 +0.040 +0.055 
-0.199 +0.197 -0.029 -0.027 -0.112 +0.949 +0.047 -0.045 +0.018 -0.034 +0.008 -0.012 
+0.197 -0.201 -0.002 +0.023 +0.010 +0.031 +0.950 -0.102 -0.027 -0.067 +0.028 +0.007 
+0.198 -0.036 -0.197 -0.032 +0.011 +0.089 +0.050 +0.949 -0.089 +0.033 +0.000 -0.003 
+0.047 -0.018 -0.198 +0.198 -0.018 -0.002 +0.022 +0.040 +0.949 +0.105 -0.069 -0.049 
+0.005 -0.029 +0.198 -0.199 -0.065 +0.029 +0.062 -0.009 -0.025 +0.949 +0.052 -0.095 
-0.013 -0.199 -0.059 +0.202 -0.026 +0.035 -0.073 -0.008 +0.014 +0.004 +0.950 +0.051 
+0.047 +0.200 -0.013 -0.196 -0.059 -0.035 +0.038 -0.000 +0.083 +0.058 +0.034 +0.949*/ 
//===================================================================================
   unitary(0,3)=alpha4;
   unitary(1,2)=alpha4;
   unitary(2,1)=alpha4;
   unitary(3,0)=alpha4;

   unitary(0,4)=+alpha1;
   unitary(2,4)=-alpha1;
   
   unitary(0,5)=+alpha1;
   unitary(1,5)=-alpha1;
   
   unitary(0,6)=-alpha1;
   unitary(1,6)=+alpha1;
   
   unitary(0,7)=-alpha1;
   unitary(2,7)=+alpha1;
   
   unitary(2,8)=+alpha1;
   unitary(3,8)=-alpha1;
   
   unitary(2,9)=-alpha1;
   unitary(3,9)=+alpha1;
   
   unitary(1,10)=+alpha1;
   unitary(3,10)=-alpha1;
   
   unitary(1,11)=-alpha1;
   unitary(3,11)=+alpha1;

   unitary(1,0)=alpha2;
   unitary(2,0)=alpha2;
   
   unitary(0,1)=alpha2;
   unitary(3,1)=alpha2;
   
   unitary(0,2)=alpha2;
   unitary(3,2)=alpha2;

   unitary(1,3)=alpha2;
   unitary(2,3)=alpha2;
   
   unitary(0,8)=+alpha3;
   unitary(0,9)=-alpha3;
   unitary(0,10)=+alpha3;
   unitary(0,11)=-alpha3;
  
   unitary(1,4)=+alpha3;
   unitary(1,7)=-alpha3;
   unitary(1,8)=-alpha3;
   unitary(1,9)=+alpha3;
   
   unitary(2,5)=+alpha3;
   unitary(2,6)=-alpha3;
   unitary(2,10)=-alpha3;
   unitary(2,11)=+alpha3;
  
   unitary(3,4)=-alpha3;
   unitary(3,5)=-alpha3;
   unitary(3,6)=+alpha3;
   unitary(3,7)=+alpha3;
   
/*Unitary
+0.913 +0.081 +0.006 +0.003 +0.200 +0.199 -0.201 -0.200 -0.028 +0.011 +0.012 -0.032 
+0.006 +0.913 -0.036 +0.088 -0.015 -0.200 +0.197 +0.035 -0.013 +0.015 +0.202 -0.198 
+0.083 +0.018 +0.913 +0.008 -0.200 +0.016 +0.001 +0.199 +0.199 -0.200 +0.055 -0.012 
-0.006 +0.006 +0.086 +0.913 -0.012 +0.041 -0.006 +0.021 -0.200 +0.198 -0.200 +0.202 
-0.199 +0.032 +0.198 -0.006 +0.949 +0.071 +0.046 +0.069 +0.070 +0.025 +0.040 +0.055 
-0.199 +0.197 -0.029 -0.027 -0.112 +0.949 +0.047 -0.045 +0.018 -0.034 +0.008 -0.012 
+0.197 -0.201 -0.002 +0.023 +0.010 +0.031 +0.950 -0.102 -0.027 -0.067 +0.028 +0.007 
+0.198 -0.036 -0.197 -0.032 +0.011 +0.089 +0.050 +0.949 -0.089 +0.033 +0.000 -0.003 
+0.047 -0.018 -0.198 +0.198 -0.018 -0.002 +0.022 +0.040 +0.949 +0.105 -0.069 -0.049 
+0.005 -0.029 +0.198 -0.199 -0.065 +0.029 +0.062 -0.009 -0.025 +0.949 +0.052 -0.095 
-0.013 -0.199 -0.059 +0.202 -0.026 +0.035 -0.073 -0.008 +0.014 +0.004 +0.950 +0.051 
+0.047 +0.200 -0.013 -0.196 -0.059 -0.035 +0.038 -0.000 +0.083 +0.058 +0.034 +0.949*/ 

 
   for (int state=0;state<nstates;state++) rotate_one_rdm(in_one_rdm[state], unitary, out_one_rdm[state]);
   double trace=0.0;
   for (int i=0;i<12;i++) trace=trace+in_one_rdm[0](i,i);
   trace=trace/4.0;
   double value=0.0;
   for (int state=0;state<nstates;state++)
   {
	   for (int i=0;i<4;i++) value=value+pow((out_one_rdm[state](i,i)-trace),2.0);
	   //value=value+pow((out_one_rdm[state](0,3)-0.0),2.0);.... This is silly
	   //value=value+pow((out_one_rdm[state](1,2)-0.0),2.0);
   }
   
   RMatrix eye;
   real_matrix_multiply_abt(unitary,unitary,eye);
   int eyecheck=4;
   for (int i=0;i<eyecheck;i++) value=value+pow((eye(i,i)-1.0),2.0);
   for (int i=0;i<eyecheck;i++) 
   {
	for (int j=0;j<eyecheck;j++)
	{
		if (i!=j) value=value+pow((eye(i,j)-0.0),2.0);
	}
   }
   //cout<<"Cost ="<<value<<endl;
   if (print_on)
   {
   	   cout<<"Cost ="<<value<<endl;
	   cout<<"Unitary"<<endl;
	   print_real_mat(unitary);
	   cout<<"Eye"<<endl;
	   print_real_mat(eye);
	   for (int state=0;state<nstates;state++) 
	   {
		cout<<"Out 1-RDM for state number "<<state<<endl;
		print_real_mat(out_one_rdm[state]);
	   }
   }
   return value;     
}

/////////////////////////////////////////////////////////////////////
void rotate_one_rdm(RMatrix &in_one_rdm,
		    RMatrix &unitary,
                    RMatrix &out_one_rdm)
{
   int nsites=in_one_rdm.NRows(); 
   out_one_rdm.resize(nsites,nsites);
   for (int i=-0;i<nsites*nsites;i++) out_one_rdm[i]=0.0;   
   RMatrix tmp;
   real_matrix_multiply(unitary,in_one_rdm,tmp); 
   real_matrix_multiply_abt(tmp,unitary,out_one_rdm);
 
 //  // Make the out one rdm as close to 0.5 on diagonals 
 //  for (int i=0;i<nsites;i++)
 //  {
 //       for (int k=0;k<nsites;k++)
 //       {
 //       	for(int j=0;j<nsites;j++)
 //       	{
 //       		for(int m=0;m<nsites;m++)
 //       		{
 //       			out_one_rdm(i,k)+=unitary(i,j)*unitary(k,m)*in_one_rdm(j,m);

 //       		}
 //       	}
 //       }
 //  }
}


void optimize_unitary_from_one_rdm(std::vector<RMatrix> &in_one_rdm, RMatrix &unitary, std::vector<RMatrix> &out_one_rdm)  
{
	double alpha1best=0.220;
	double alpha2best=0.044;
	double alpha3best=0.020;
	double alpha4best=0.018;
        /*double mincost=1000; 
	for (int i1=0;i1<=125;i1++)
	{
		double alpha1=double(i1)/500.0;
		for (int i2=0;i2<=50;i2++)
		{
			double alpha2=double(i2)/500.0;
			for (int i3=-50;i3<=50;i3++)
			{
				double alpha3=double(i3)/500.0;
				for (int i4=-50;i4<=50;i4++)
				{
					  double alpha4=double(i4)/500.0;
					  double c=cost(in_one_rdm, alpha1, alpha2, alpha3, alpha4, false);
					  if (c<mincost) 
					  {
						mincost=c;
						alpha1best=alpha1;
						alpha2best=alpha2;
						alpha3best=alpha3;
						alpha4best=alpha4;
					  }
				}
			}
		}
	}*/
	double c=cost(in_one_rdm, alpha1best, alpha2best, alpha3best, alpha4best, true);
}
