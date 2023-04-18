#include <stdio.h>
#include <math.h>

// CONSTANTS
double e = 0.0167;            // ECCENTRICITY
double eps = 0.4089889;       // OBLICUITY OF ORBIT
double Lp = 4.943063;         // ECLIPTIC LONGITUDE OF PERIGEE

double TOL = 1e-10;
double year = 365.2425;

char buffer[50];
int cal[12]={31,28,31,30,31,30,
             31,31,30,31,30,31};
char* monthName[12]={"jan","feb","mar",
                     "apr","may","jun",
                     "jul","aug","sep",
                     "oct","nov","dec"};


double kepler(double E, double M) {
  return E-e*sin(E)-M;
}
double Dkepler(double E, double M) {
  return 1-e*cos(E);
}

double mod(double op1, double op2) {
  while(op1<0) op1+=op2;
  while(op1>=op2) op1-=op2;
  return op1;
}

void tic(int *day,int *month){
  *day=*day+1;
  if(*month == 12 && *day == 32) {
    *day=*month=1;
    return;
  }
  if(*day>cal[*month-1]){
    *day=1;
    *month=*month+1;
    return;
  }
}

char* dateToString(int n){
  int month=1;
  int day=1;
  int i;
  for(i=0;i<n;++i){
    tic(&day,&month);
  }
  sprintf(buffer,"%s, %2i%s",
        monthName[month-1],
        day,
        (day==1)?"st":((day==2)?"nd":"th"));
  return buffer;
}

double newton(double (*f)(double,double), 
  double (*df)(double,double),
  double E, double M) {
  
  double newE=E;
  do {
    E=newE;
    newE = E-((*f)(E,M))/((*df)(E,M));
  } while (fabs(E-newE)>TOL);
  return newE;
}

int main(int argc, char const *argv[]) {
  int day;            // DAY OF YEAR
  double M;           // MEAN ANOMALY OF REAL SUN
  double E;           // ECCENTRIC ANOMALY OF REAL SUN
  double phi;         // REAL ANOMALY OF REAL SUN
  double L;           // ECLIPTIC LONGITUE OF REAL SUN
  double alpha;       // RIGHT ASCENSION OF REAL SUN
  double alphaM;      // RIGHT ASCENSION OF MEAN SUN
  double eot;         // ÆQUATION OF TIME

  int month;

  // JAN 1st
  day=74;
  M = mod(((day-4.0)/year*2*M_PI),2*M_PI);
  E = M;
  E = mod(newton(&kepler,&Dkepler,E,M),2*M_PI);
  phi = mod(2*atan(sqrt((1+e)/(1-e))*tan(E/2)),2*M_PI);
      // RECHECK ANGLE? I DON'T THINK SO
  L=mod(Lp+phi,2*M_PI);
  alpha = atan(cos(eps)*tan(L));
  while(alpha-L < 1.0) alpha+=M_PI;
  while(alpha-L > 1.0) alpha-=M_PI;


  fprintf(stderr,"M = %e\n", M);
  fprintf(stderr,"E = %e\n", E);
  fprintf(stderr,"phi = %e\n",phi);
  fprintf(stderr,"L = %e\n",L);
  fprintf(stderr,"alpha = %e\n\n",alpha);



  alphaM = mod(Lp+M,2*M_PI);
  eot = alphaM - alpha;


  fprintf(stderr,"alphaM = %e\n\n",alphaM);
  fprintf(stderr,"eot (min) = %e\n\n",eot/(M_PI)*12*60);


  //for(month=0,day=0;day<365;day+=cal[month++]) {
  for(day=0;day<365;day++) {
    M = mod(((day-4.0)/year*2*M_PI),2*M_PI);
    E = M;
    E = mod(newton(&kepler,&Dkepler,E,M),2*M_PI);
    phi = mod(2*atan(sqrt((1+e)/(1-e))*tan(E/2)),2*M_PI);
    L=mod(Lp+phi,2*M_PI);
    alpha = atan(cos(eps)*tan(L));
    while(alpha-L < 1.0) alpha+=M_PI;
    while(alpha-L > 1.0) alpha-=M_PI;
    
    alphaM = mod(Lp+M,2*M_PI);
    eot = alphaM - alpha;
    while(eot>M_PI/2) eot-=M_PI;
    while(eot<-M_PI/2) eot+=M_PI;
    //printf("\"%s\",\t%+f\n",dateToString(day), eot/(M_PI)*12*60);
    printf("%3i %+f\n",day, eot/(M_PI)*12*60);
  }

}