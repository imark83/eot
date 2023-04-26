#include <stdio.h>
#include <math.h>

// CONSTANTS
double e = 0.01669;           // ECCENTRICITY
double eps = 0.4089889;       // OBLICUITY OF ORBIT
// double Lp = 4.943063;         // ECLIPTIC LONGITUDE OF PERIGEE
double Lp = 4.975;            // ECLIPTIC LONGITUDE OF PERIGEE
double latZgz = 0.726958722276505; // GEOGRAPHICAL LATITUDE OF ZGZ

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

int main(int argc, char const *argv[]){
  int day;            // DAY OF YEAR
  double M;           // MEAN ANOMALY OF SUN
  double E;           // ECCENTRIC ANOMALY OF SUN
  double phi;         // TRUE ANOMALY OF SUN
  double L;           // ECLIPTIC LONGITUE OF SUN
  double delta;       // DECLINATION OF SUN
  double halfDay;     // LENGTH OF HALF OF THE DAY

  double hours, minutes;
  double aux;

  int month;


  day=79;
  M = mod(((day-4.678472)/year*2*M_PI),2*M_PI);
  E = M;
  E = mod(newton(&kepler,&Dkepler,E,M),2*M_PI);
  phi = mod(2*atan(sqrt((1+e)/(1-e))*tan(E/2)),2*M_PI);
      // RECHECK ANGLE? I DON'T THINK SO
  L=mod(Lp+phi,2*M_PI);
  delta = asin(sin(eps)*sin(L));
  halfDay = acos(-tan(latZgz)*tan(delta));
  aux = modf(2*halfDay/M_PI*12,&hours);
  modf(aux*60,&minutes);

  printf("L = %f\n",L);
  printf("delta = %f\n",delta);
  printf("halfDay = %f\n",halfDay);

  printf("%s -> %i hours, %i minutes\n",
          dateToString(day),(int) hours,(int) minutes);



  for(day=0;day<365;day++) {
    M = mod(((1+day-4.678472)/year*2*M_PI),2*M_PI);
    E = M;
    E = mod(newton(&kepler,&Dkepler,E,M),2*M_PI);
    phi = mod(2*atan(sqrt((1+e)/(1-e))*tan(E/2)),2*M_PI);
        // RECHECK ANGLE? I DON'T THINK SO
    L=mod(Lp+phi,2*M_PI);
    delta = asin(sin(eps)*sin(L));
    halfDay = acos(-tan(latZgz)*tan(delta));
    aux = modf(2*halfDay/M_PI*12,&hours);
    modf(aux*60,&minutes);

    delta*=180/M_PI;
    L*=180/M_PI;
    phi*=180/M_PI;
    M*=180/M_PI;
    printf("%s -> %2i hours, %2i minutes\n",
           dateToString(day),(int) hours,(int) minutes);

    // printf("%s\tphi=%6.2f\tM=%6.2f\tdelta=%6.2f\tL=%6.2f\n",
    //         dateToString(day),phi,M,delta,L);


  }

  return 0;
}