#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <memory.h>
#include <fstream>
#include <iostream>

#define DUMP	1
#define MAX_FFA	1000
#define MAX_D	1000

using namespace std;

int D = 30;					// Sorunun boyutu
int n = 20;						// Ateþböceði sayýsý
int MaxGeneration = 1000;			// Yineleme(iterasyon) sayýsý
int NumEval;					// Deðerlendirme sayýsý
int Index[MAX_FFA];				// Fitness deðerlerine göre ateþböcekleri çeþitleri

double ffa[MAX_FFA][MAX_D];		// firefly sayisi
double ffa_tmp[MAX_FFA][MAX_D]; // ara populasyon
double f[MAX_FFA];				//uygunluk
double I[MAX_FFA];				// Iþýk þiddeti
double nbest[MAX_FFA];          // en iyi cozum
double lb[MAX_D];				// Üst sýnýr
double ub[MAX_D];				// Alt sýnýr

double alpha = 0.5;				// alpha parametresi //0.5
double betamin = 0.2;           // beta parametresi //0.2
double gama = 1.0;				// gamma parametresi //1.0

double fbest;					// En iyi amaç fonksiyonu

typedef double (*FunctionCallback)(double sol[MAX_D]);

/*benchmark fonksiyonlarý */
double sphere(double sol[MAX_D]);
double sumSquares(double sol[MAX_D]);

/*Kendi amaç fonksiyonunu yazýn */
FunctionCallback function = &sphere; //////////////////

// Isteðe baðlý olarak yeni alfa deðerini tekrar hesapla
double alpha_new(double alpha, int NGen) ///////////////////
{
	double delta;			// delta parameter
	delta = 1.0-pow((pow(10.0, -4.0)/0.9), 1.0/(double) NGen);
	return (1-delta)*alpha;
}

// Ateþböceði popülasyonunu baþlat
void init_ffa() 
{
	int i, j;
	double r;

	// Üst ve alt sýnýrlarý baþlat
	for (i=0;i<D;i++)
	{
		lb[i] = -100.0;
		ub[i] = 100.0;
	}

	for (i=0;i<n;i++)
	{
		for (j=0;j<D;j++)
		{
			r = (   (double)rand() / ((double)(RAND_MAX)+(double)(1)) );
			ffa[i][j]=r*(ub[j]-lb[j])+lb[j];
		}
		f[i] = 1.0;			// Çekiciliði baþlatýn
		I[i] = f[i];		//Iþýk þiddetlerini atadýk
	}
}

// Kabarcýk sýralama uygulamasý
void sort_ffa() 
{
	int i, j;

	// Endekslerin baþlatýlmasý
	for(i=0;i<n;i++)
	{
		Index[i] = i;
	}		

	// Kabarcýk sýralama
	for(i=0;i<n-1;i++)
	{
		for(j=i+1;j<n;j++)
		{
			if(I[i] > I[j])
			{
				double z = I[i];	// Deðiþim cazibesi
				I[i] = I[j];
				I[j] = z;
				z = f[i];			// Deðiþim alýþkanlýðý
				f[i] = f[j];
				f[j] = z;
				int k = Index[i];	// Döviz endeksleri
				Index[i] = Index[j];
				Index[j] = k;
			}
		}
	}
}

// Eski popülasyonu yeni Ýndeks deðerlerine göre deðiþtirin
void replace_ffa()
{
	int i, j;

	// Orijinal nüfusu geçici bölgeye kopyala
	for(i=0;i<n;i++)
	{
		for(j=0;j<D;j++)
		{
			ffa_tmp[i][j] = ffa[i][j];
		}
	}

	// EA anlamýnda nesil seçimi
	for(i=0;i<n;i++)
	{
		for(j=0;j<D;j++)
		{
			ffa[i][j] = ffa_tmp[Index[i]][j];
		}
	}
}

void findlimits(int k) //lb ve ub aralýðýndakileri seçme
{
	int i;

	for(i=0;i<D;i++)
	{
		if(ffa[k][i] < lb[i])
			ffa[k][i] = lb[i];
		if(ffa[k][i] > ub[i])
			ffa[k][i] = ub[i];
	}
}

void move_ffa() 
{
	int i, j, k;
	double scale;
	double r, beta;

	for(i=0;i<n;i++)
	{
		scale = abs(ub[i]-lb[i]);
		for(j=0;j<n;j++)
		{
			r = 0.0;
			for(k=0;k<D;k++)
			{
				r += (ffa[i][k]-ffa[j][k])*(ffa[i][k]-ffa[j][k]);
			}
			r = sqrt(r);
			if(I[i] > I[j])	// Daha parlak ve daha cazip
			{
				double beta0 = 1.0;
				beta = (beta0-betamin)*exp(-gama*pow(r, 2.0))+betamin;
				for(k=0;k<D;k++)
				{
					r = (   (double)rand() / ((double)(RAND_MAX)+(double)(1)) );
					double tmpf = alpha*(r-0.5)*scale;
					ffa[i][k] = ffa[i][k]*(1.0-beta)+ffa_tmp[j][k]*beta+tmpf;
				}
			}
		}
		findlimits(i);
	}
}

void dump_ffa(int gen) 
{
	cout << "Iterasyon: " << gen << " En iyi: " << fbest << endl;
	
	
}



int main(int argc, char* argv[]) 
{
    ofstream yaz;
      
    yaz.open("veriler.txt");

	int i;
    int t = 1;		// Nesil sayaç
	 
	// firefly algoritma optimizasyon döngüsü
	// rastgele üreticinin baþlangýç noktasýný belirler
	srand(1);

	//N ateþ böreklerinin ilk yerlerini üretmek
	init_ffa();
	
	#ifdef DUMP
	
	dump_ffa(t);
	
	yaz << "Iterasyon: " << t;
	yaz << "		";
	yaz << "Best: " <<fbest;
	yaz << "\n";
	
	#endif

	while(t <= MaxGeneration)
	{
		// Alfa küçültme bu satýrý isteðe baðlýdýr
		alpha = alpha_new(alpha, MaxGeneration);

		// Yeni çözümleri deðerlendirmek
		for(i=0;i<n;i++)
		{
            f[i] = function(ffa[i]);         // Çözümün uygunluðunu saðlamak. Probleme ateþböceklerini yolluyoruz.
			I[i] = f[i];					// Çekiciliði baþlatýn
		}

		// Iþýk þiddetine göre ateþböceklerini sýralamak
		sort_ffa();
		// Eski nüfusun yerini al
		replace_ffa();

		// En iyi ürünü bul
		for(i=0;i<D;i++)
	    {
			nbest[i] = ffa[0][i];
	    }			
		fbest = I[0];

		// Tüm ateþböceklerini daha iyi yerlere taþýyýn
		move_ffa();
		#ifdef DUMP
		dump_ffa(t);
		#endif
		t++;
		
		yaz << "Iterasyon: " << t;
		yaz << "		";
		yaz << "Best: " <<fbest;
		yaz << "\n";
   
	}
	
	yaz.close();
	cout << "Optimizasyon sonu: En iyisi: " << fbest << endl;

	return 0;
}

double sphere(double* sol)  //Sphere Fonksiyonu
{
	int j;
	double top = 0;
	for (j = 0; j < D; j++) 
	{
		top = top + (sol[j] * sol[j]);
	}
	return top;
}

double sumSquares(double* sol) //Sum Squares Fonksiyonu
{
	int j;
	double top = 0;
	for (j = 0; j < D; j++) 
	{
		top = top + (j * (sol[j] * sol[j]));
	}
	return top;
}
