// Ising2D T3.cpp
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#define N	26	// Trebuie pus cu 1 mai mult
#define JJ	1
#define T	3
#define Kb	1
#define NrPasiMonterCarlo 10001			// Trebuie pus cu 1 mai mult
#define _CRT_SECURE_NO_WARNINGS

void numarare(int a[N][N])
{
	int n = N;
	for (int i = 1; i <= n; i++)
	{
		int nr = 0;
		for (int j = 1; j <= n; j++)
			if (a[i][j] == 1)
				nr++;
			else
			{
				if (nr != 0)
					printf("linia, numar+ %d,%d\n", i, nr);
				nr = 0;
			}
	}
}

int corelatie(int spin[N][N])
{
	FILE* fptr;
	fptr = fopen("C.txt", "w");
	int C = 0;
	int Suma = 0;
	for (int x = 1; x < N - 1; x++)
	{
		for (int y = 1; y < N - 1; y++)
		{
			Suma = spin[x - 1][y] + spin[x + 1][y] + spin[x][y - 1] + spin[x][y + 1];
			C += Suma * spin[x][y];
			// Suma = 0;
		}
	}
	printf("%d", C);
	fprintf(fptr, "%d", C);
	fclose(fptr);
	return C;
}

int corelatie_omogena(int spin[7][7])
{
	FILE* fptr;
	fptr = fopen("C.omogen.txt", "w");
	int C = 0;
	int Suma = 0;
	for (int x = 1; x < 7 - 1; x++)
	{
		for (int y = 1; y < 7 - 1; y++)
		{
			Suma = spin[x - 1][y] + spin[x + 1][y] + spin[x][y - 1] + spin[x][y + 1];
			C += Suma * spin[x][y];
			// Suma = 0;
		}
	}
	printf("%d", C);
	fprintf(fptr, "%d", C);
	fclose(fptr);
	return C;
}

void matrice_omogena(int a[7][7], int o)
{
	for (int i = 0; i < o; i++)// Atribuim tuturor spinilor 1
	{
		for (int j = 0; j < o; j++)
		{
			a[i][j] = 1;
		}
	}
	for (int j = 0; j < o; j++)// Egalam marginile cu zero
	{
		a[0][j] = 0;
		a[o - 1][j] = 0;
		a[j][0] = 0;
		a[j][o - 1] = 0;
	}
	for (int i = 1; i < o - 1; i++)
		for (int j = 1; j < o - 1; j++)
			if (i < 3)
				a[i][j] = 1;
			else
				a[i][j] = -1;
}

int get_energy(int b[N][N])
{
	int a[N][N];

	for (int i = 0; i < N; i++)
		for (int j = 0; j < N; j++)
			a[i][j] = b[i][j];

	int Sv = 0;
	int dE = 0;

	for (int i = 1; i < N - 2; i++)
		for (int j = 1; j < N - 2; j++) {
			//			if (i >= 1 && i <= N - 2 && j >= 1 && j <= N - 2)  // Fara margini
			//			{
			Sv = a[i - 1][j] + a[i + 1][j] + a[i][j - 1] + a[i][j + 1];	// Calculam diferenta de energie
			dE += 2 * JJ * Sv * a[i][j];
			//printf("%d\n", dE);
//			}
		}

	return dE;
}

int tryGetMinEnergy(int a[N][N])
{
	int x = 0;
	int y = 0;
	int minEnergy = get_energy(a);
	FILE* fptr = fopen("EnergiiMinime.txt", "w");
	for (int i = 1; i < N - 1; i++)
	{
		for (int j = 1; j < N - 1; j++)
		{
			for (int direction = 0; direction <= 3; direction++)
			{

				if (direction == 0)
				{
					x = i;
					y = j + 1;
				}
				if (direction == 1)
				{
					x = i + 1;
					y = j;
				}
				if (direction == 2)
				{
					x = i;
					y = j - 1;
				}
				if (direction == 3)
				{
					x = i - 1;
					y = j;
				}

				if (x >= 1 && x <= N - 2 && y >= 1 && y <= N - 2)
				{
					if (a[i][j] != a[x][y])
					{
						int unchangedEnergy = get_energy(a);

						a[i][j] *= -1;
						a[x][y] *= -1;

						int changedEnergy = get_energy(a);

						//printf("\nEnergii: %d versus %d\n", unchangedEnergy, changedEnergy);
						//printf("\nPentru valoarea a[%d][%d]=%d verificam vecinul a[%d][%d]=%d, energia minima=%d, matricea:\n", i, j, a[i][j]*-1, x, y, a[x][y]*-1, changedEnergy);
						//for (int i = 0; i < N; i++)
						//{
						//	for (int j = 0; j < N; j++)
						//	{
						//		printf("%2d ", a[i][j]);
						//	}
						//	printf("\n");
						//}

						if (changedEnergy < unchangedEnergy)
						{
							minEnergy = changedEnergy;
							i = 1;
							j = 1;
							break;
						}
						else
						{
							a[i][j] *= -1;
							a[x][y] *= -1;
						}
					}
				}
			}
		}

	/*	fprintf(fptr, "%d\n", minEnergy);*/

	}
	fclose(fptr);

	return minEnergy;
}

int spin[N][N];

void generareGlauber() 
{
	int i, j;
	int Sv = 0;	// Suma spinilor vecini - Nord Sud Est Vest
	int M[NrPasiMonterCarlo] = { 0 };// Magnetizatia
	double P;
	double r;
	int dE = 0;
	FILE* fptr;

	for (i = 0; i < N; i++)// Atribuim tuturor spinilor 1
	{
		for (j = 0; j < N; j++)
		{
			spin[i][j] = 1;
		}
	}


	for (j = 0; j < N; j++)// Egalam marginile cu zero
	{
		spin[0][j] = 0;
		spin[N - 1][j] = 0;
		spin[j][0] = 0;
		spin[j][N - 1] = 0;
	}


	int k = 0;
	srand(time(NULL)); // initializezi cu timpu de pe procesor
	for (int pasiMC = 0; pasiMC < NrPasiMonterCarlo; pasiMC++)// Pasi Monte Carlo
	{
		for (int contor1 = 0; contor1 < N * N; contor1++)
		{
			int x = int(rand() / (RAND_MAX + 1.0) * N);  // Alegem spinii carora vrem sa le schimbam semnul
			int y = int(rand() / (RAND_MAX + 1.0) * N);

			if (x >= 1 && x <= N - 2 && y >= 1 && y <= N - 2)  // Fara margini
			{
				Sv = spin[x - 1][y] + spin[x + 1][y] + spin[x][y - 1] + spin[x][y + 1];	// Calculam diferenta de energie
				dE = 2 * JJ * Sv * spin[x][y];
				P = exp(-dE / T) / (1 + exp(-dE / T));// Probabilitatea ca spinii sa se schimbe	Beta =
				r = (rand() / (RAND_MAX + 1.0));
				if (r < P)
				{
					spin[x][y] = -spin[x][y];
				}
			}

		}
	}
	//scrie matricea spin in fisier
	fptr = fopen("Glauber Patrat T0.3.txt", "w");
	for (i = 0; i < N; i++)
	{
		for (j = 0; j < N; j++)
		{
			fprintf(fptr, "\n");
			fprintf(fptr, "%2d\t%2d\t%2d\n ", i, j, spin[i][j]);
		}
	}
	fclose(fptr);
	//afiseara matricea spin pe ecran
	printf("\n");
	/*for (i = 0; i < N; i++)
	{
		for (j = 0; j < N; j++)
		{
			printf("%2d ", spin[i][j]);
		}
		printf("\n");
	}*/
}

void generareMetropolis() {
	// Matricea
	int i, j;
	int Sv = 0;	// Suma spinilor vecini - Nord Sud Est Vest
	int M[NrPasiMonterCarlo] = { 0 };// Magnetizatia
	double P;
	double r;
	int dE = 0;
	FILE* fptr;

	for (i = 0; i < N; i++)// Atribuim tuturor spinilor 1
	{
		for (j = 0; j < N; j++)
		{
			spin[i][j] = 1;
		}
	}


	for (j = 0; j < N; j++)// Egalam marginile cu zero
	{
		spin[0][j] = 0;
		spin[N - 1][j] = 0;
		spin[j][0] = 0;
		spin[j][N - 1] = 0;
	}


	int k = 0;
	srand(time(NULL));
	for (int pasiMC = 0; pasiMC < NrPasiMonterCarlo; pasiMC++)// Pasi Monte Carlo
	{
		for (int contor1 = 0; contor1 < N * N; contor1++)
		{
			int x = int(rand() / (RAND_MAX + 1.0) * N);  // Alegem spinii carora vrem sa le schimbam semnul
			int y = int(rand() / (RAND_MAX + 1.0) * N);

			if (x >= 1 && x <= N - 2 && y >= 1 && y <= N - 2)  // Fara margini
			{
				Sv = spin[x - 1][y] + spin[x + 1][y] + spin[x][y - 1] + spin[x][y + 1];	// Calculam diferenta de energie
				dE = 2 * JJ * Sv * spin[x][y];
				P = exp(-dE / T);// Probabilitatea ca spinii sa se schimbe	Beta =
				r = (rand() / (RAND_MAX + 1.0));
				if (r < P)
				{
					spin[x][y] = -spin[x][y];
				}
			}

		}
	}
	//scrie matricea spin in fisier
	fptr = fopen("Metropolis Patrat T0.3.txt", "w");
	for (i = 0; i < N; i++)
	{
		for (j = 0; j < N; j++)
		{
			fprintf(fptr, "\n");
			fprintf(fptr, "%2d\t%2d\t%2d\n ", i, j, spin[i][j]);
		}
	}
	fclose(fptr);
	//afiseara matricea spin pe ecran
	printf("\n");
	/*for (i = 0; i < N; i++)
	{
		for (j = 0; j < N; j++)
		{
			printf("%2d ", spin[i][j]);
		}
		printf("\n");
	}*/
}

int main()
{

	generareMetropolis();
	for (int i = 0; i < N; i++)
	/*{
		for (int j = 0; j < N; j++)
		{
			printf("%2d ", spin[i][j]);
		}
		printf("\n");
	}*/

	//numara spinii de acelasi semn aflati unu langa celalalt pe aceeasi linie
	/*numarare(spin);*/

	//calculeaza factorul de corelatie
	/*corelatie(spin);*/

	//initializam matricea omogena
	/*int a[7][7];
	matrice_omogena(a, 7);*/

	//afisam matricea omogena
	//printf("\n");
	//for (i = 0; i < 7; i++)
	//{
	//	for (j = 0; j < 7; j++)
	//	{
	//		printf("%2d ", a[i][j]);
	//	}
	//	printf("\n");
	//}
	// afisam pe ecran corelatia matricei omogena
	//corelatie_omogena(a);

	//scriem in fisier matricea omogena
	/*fptr = fopen("Omogena.txt", "w");
	for (i = 0; i < N; i++)
	{
		for (j = 0; j < N; j++)
		{
			fprintf(fptr, "\n");
			fprintf(fptr, "%2d\t%2d\t%2d\n ", i, j, a[i][j]);
		}
	}
*/

	printf("\n energie: %d\n", get_energy(spin));
	printf("\n energie minima: %d", tryGetMinEnergy(spin));

	/*printf("\nMatricea finala cu energia cea mai mica va fi:\n");*/
	/*for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			printf("%2d ", spin[i][j]);
		}
		printf("\n");
	}*/

	generareGlauber();

	/*for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			printf("%2d ", spin[i][j]);
		}
		printf("\n");
	}*/

	printf("\n energie: %d\n", get_energy(spin));
	printf("\n energie minima: %d", tryGetMinEnergy(spin));

	return 0;
}
