#include "mbed.h"
#include <complex>
#include "math.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#define MAX 2048

#define M_PI 3.1415926535897932384

using namespace std;

Serial pc(USBTX,USBRX);


int log2(int N)    /*function to calculate the log2(.) of int numbers*/
{
  int k = N, i = 0;
  while(k) {
    k >>= 1;
    i++;
  }
  return i - 1;
}

int check(int n)    //checking if the number of element is a power of 2
{
  return n > 0 && (n & (n - 1)) == 0;
}

int reverse(int N, int n)    //calculating revers number
{
  int j, p = 0;
  for(j = 1; j <= log2(N); j++) {
    if(n & (1 << (log2(N) - j)))
      p |= 1 << (j - 1);
  }
  return p;
}

void ordina(complex<float>* f1, int N) //using the reverse order in the array
{
  complex<float> f2[MAX];
  for(int i = 0; i < N; i++)
    f2[i] = f1[reverse(N, i)];
  for(int j = 0; j < N; j++)
    f1[j] = f2[j];
}

void transform(complex<float>* f, int N) //
{
  ordina(f, N);    //first: reverse order
  complex<float> *W;
  W = (complex<float> *)malloc(N / 2 * sizeof(complex<float>));
  W[1] = polar(1., -2. * M_PI / N);
  W[0] = 1;
  for(int i = 2; i < N / 2; i++)
    W[i] = pow(W[1], i);
  int n = 1;
  int a = N / 2;
  for(int j = 0; j < log2(N); j++) {
    for(int i = 0; i < N; i++) {
      if(!(i & n)) {
        complex<float> temp = f[i];
        complex<float> Temp = W[(i * a) % (n * a)] * f[i + n];
        f[i] = temp + Temp;
        f[i + n] = temp - Temp;
      }
    }
    n *= 2;
    a = a / 2;
  }
}

void FFT(complex<float>* f, int N, float d)
{
  transform(f, N);
  //for(int i = 0; i < N; i++)
    //f[i] *= d; //multiplying by step
}

int main()
{
  srand(time(NULL));
  int n = MAX;
  /*do {
    pc.printf("Specifiy array dimension (MUST be a power of 2)\n\r");
    pc.scanf("%d",&n);
  } while(!check(n)); */
  float d = 1;
  //pc.printf("Specify the sampling step\n\r");
  //pc.scanf("%lf",&d);
  
  complex<float> vec[MAX];
  int harmonics = 0;
  float Fs = 5000;
  float T = 1/Fs;
  float t;  
  complex<float> vecc[MAX];
  float S[MAX];
  float positiveband[(MAX/2)+1];
  float frequencies[(MAX/2)];
  
  /*Riempimento del vettore di misure*/
  for(int i = 0; i < n; i++) {
    t = i*T;
    vec[i] = complex<float>(10*sin(2*M_PI*50*t)/*+5*cos(2*M_PI*200*t)+sin(2*M_PI*300*t)*/,0);
    //wait_ms(10); //prende un campione ogni 10 ms - 100 Hz    
  }
   
  clock_t begin = clock();
  
  /*Fa la FFT del vettore campionato*/
  FFT(vec, n, d);    
  
  /*Fa la normalizzazione dei campioni dividendoli per la lunghezza del vettore*/
  for(int j = 0; j < n; j++){
      vecc[j] = complex<float>(vec[j].real()/MAX,vec[j].imag()/MAX); 
  }  
  
  /*Fa il valore assoluto dei campioni*/
  for(int j = 0; j < n; j++){
      S[j] = abs(vecc[j]); 
  }  
  
  /*Prende solo la parte positiva dello spettro*/
  for(int j = 1; j < (MAX/2)+1; j++){
      positiveband[j] = S[j]; 
  }
  
  /*Raddoppia la parte positiva dello spettro*/
  for(int j = 2; j < (MAX/2)+1; j++){
      positiveband[j] = 2*positiveband[j]; 
  }
  
  int location;
  float maximum;
  
  maximum = positiveband[0];
  for (int j = 1; j < (MAX/2)+1; j++){
    if (positiveband[j] > maximum){
        maximum  = positiveband[j];
        location = j;
    }
  }
  
  for(int j = 0; j < (MAX/2); j++){
      frequencies[j] = (Fs*j)/MAX; 
  }
  
  int i;
  for(i = 1; i <= 50; i++){
       pc.printf("%d harmonics, amplitude: %f, frequency: %f\n\r",i,positiveband[i*location],frequencies[i*location]);
  }  
  
  float total = 0.0; 
  
  for(i = 2; i <= 50; i++){
       //pc.printf("%d harmonics, amplitude: %f, frequency: %f\n\r",i,positiveband[i*location],frequencies[i*location]);
       total += positiveband[i*location]*positiveband[i*location];
  }  
  
    float thd = (sqrt(total)/positiveband[location])*100;
    
    pc.printf("THD: %f\n\r",thd);
    
   
  //Stampa le ampiezze 
  /*for(int j = 0; j < (MAX/2)+1; j++){
      pc.printf("positiveband[%d]: %f\n\r",j,positiveband[j]);      
  } */

    
  //Stampa le frequenze: deve arrivare fino a 2500 HZ perché 50*f0
  /*for(int j = 0; j < (MAX/2); j++){
      pc.printf("frequencies[%d]: %f\n\r",j,frequencies[j]);       
  } */
  
   
  clock_t end = clock();
  float time_spent = (float)(end - begin) / CLOCKS_PER_SEC;
  
  pc.printf("Elapsed time (seconds): %f for %d samples \n\r",time_spent,MAX);
  
  return 0;
}


