#include<iostream>
#include<cmath>
#include<fstream>

using namespace std;

void step(double*,double*,const double, const int);
void hamiltonian(double*,double*,double&);

int main(){

const int dim=2;
double p[dim];
double q[dim];

const double e=0.6;
const double pi=3.141;
const double tEnd=20*pi;
const double t0=0;
const double dt= 0.0005;
const double N=(tEnd-t0)/dt;
double t=t0;



q[0]=1-e;
q[1]=0;
p[0]=0;
p[1]=sqrt((1+e)/(1-e));

double H;
hamiltonian(q,p,H);

ofstream out ("step.txt");
out <<t<<"\t"<<H<<"\t"<< q[0]<<"\t"<<q[1]<<endl;

for(int i=0;i<N;i++){
  step(q,p,dt,dim);
  hamiltonian(q,p,H);
  t+=dt;
  
out <<t<<"\t"<<H<<"\t"<< q[0]<<"\t"<<q[1]<<endl;
}
out.close();
return 0 ;
}

void step(double* q, double* p,const double dt,const int dim){
   double Q=q[0]*q[0]+q[1]*q[1];
  for (int i =0;i<dim;i++){
  
  p[i]-=dt*(q[i]/(pow(Q,1.5)));
  q[i]+=dt*p[i];
  }
}
  
void hamiltonian(double* q,double* p,double& H){
H=0.5*(p[0]*p[0]+p[1]*p[1])-1.0/sqrt(q[0]*q[0]+q[1]*q[1]);
}
