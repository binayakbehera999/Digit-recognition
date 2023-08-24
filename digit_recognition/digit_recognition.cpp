// digit_recognition.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <stdio.h>
#include <conio.h>
#include <stdlib.h>
#include <math.h>
#define MAX 320
#define PI 3.14
#define delta 0.001
#define epsilon 0.03
#define T_MAX 120
#define N 5
#define M 32
#define MIN 1e-30

char inputFile1[] = "English_Modified/224101014_E_0_0.txt";
char inputFile2[] = "English_Modified/224101014_E_0_10.txt";
char test[] = "English_Modified/224101014_E_0_20.txt";
char tempFile1[] = "temp/224101014_E_0_0.txt";
char tempFile2[] = "temp/224101014_E_0_10.txt";
char tempFile3[] = "temp/224101014_E_0_20.txt";
char finalFile1[] = "FinalModel/224101014_E_0_A.txt";
char finalFile2[] = "FinalModel/224101014_E_0_B.txt";

FILE *uni, *logfp;

int line = 0;
int v_count[32] = {0}, l_index = 0;
long double DC_Shift = 0, tokura_wt[12] = {1.0, 3.0, 7.0, 13.0, 19.0, 22.0, 25.0, 33.0, 42.0, 50.0, 56.0, 61.0};
long double n_amp[MAX];
long double alpha[13][13];
long double a_lpc[13];
long double r[13];
long double c[13];
long double dist_hist[1000] = {0.0};
long double centroid[32][12] = {0.0};
long double codebook[32][12] = {0.0};
long double distortion[32] = {0};
long double distance[32] = {0};
long double universeIn[9106][12];
int o[T_MAX];
int T;
long double a[N + 1][N + 1], b[N + 1][M + 1], pi[N + 1];
long double alphaf[T_MAX][N + 1];
long double beta[T_MAX][N + 1];
long double gamma[T_MAX][N + 1];
long double del[T_MAX][N + 1];
int si[T_MAX][N + 1];
long double theta[T_MAX][N + 1];
long double zeta[T_MAX][N + 1][N + 1];
long double p_star;
int q_star[T_MAX];
long double p_new = 0;
long double p_old = 0;
long double a_dash[N + 1][N + 1] = {0.0};
long double b_dash[N + 1][M + 1] = {0.0};
long double pTest[10];

// @desc Calculate DC Shift for the input file
void calculate_DCShift(char *filename)
{
	long double sum_Amp = 0, count = 0;
	FILE *fp;
	fp = fopen(filename, "r");
	if (!fp)
		perror("fp not present");

	int current;
	while (!feof(fp))
	{
		fscanf(fp, "%d", &current);
		sum_Amp = sum_Amp + abs(current);
		count++;
	}
}

// @desc Find the Max Amplitude for Normalisation
long double find_Max_Amp(long double amp[])
{
	long Max = amp[0];
	for (int i = 1; i < MAX; i++)
	{
		if (fabs(amp[i]) > Max)
			Max = fabs(amp[i]);
	}
	return Max;
}

// @desc print normalised amplitude data into output file
void normalise_Amplitude()
{

	long double normalisation_Factor;
	long double max_Amp = find_Max_Amp(n_amp);
	normalisation_Factor = 5000.0 / max_Amp;
	for (int j = 0; j < MAX; j++)
	{
		n_amp[j] = (n_amp[j] - DC_Shift) * normalisation_Factor;
	}
}

// @desc apply hamming window on the input
void hammingWindow()
{
	long double w;
	for (int j = 0; j < MAX; j++)
	{
		w = 0.54 - 0.46 * (cos((2 * PI * j) / (319.0)));
		n_amp[j] = n_amp[j] * w;
	}
}

// Calculate ai's using durbin's algorithm
int durbinAlgo()
{
	long double sum = 0;
	long double e[13], k[13];
	e[0] = r[0];
	for (int i = 1; i < 13; i++)
	{

		for (int j = 1; j <= i - 1; j++)
		{
			sum = sum + alpha[i - 1][j] * r[i - j];
		}

		k[i] = (r[i] - sum) / e[i - 1];

		alpha[i][i] = k[i];

		for (int j = 1; j <= i - 1; j++)
		{
			alpha[i][j] = alpha[i - 1][j] - (k[i] * alpha[i - 1][i - j]);
		}

		e[i] = (1.0 - (k[i] * k[i])) * e[i - 1];

		sum = 0;
	}

	for (int x = 1; x < 13; x++)
	{
		a_lpc[x] = alpha[12][x];
	}

	return 1;
}

// Calculate cepstal coefficients
void calculateCepstal()
{
	double t, sw = 0;
	for (int m = 1; m < 13; m++)
	{
		c[m] = a_lpc[m];
		for (int j = 1; j <= m - 1; j++)
		{
			t = double(j) / double(m);

			c[m] = c[m] + t * c[j] * a_lpc[m - j];
		}
	}

	// apply sine window
	for (int m = 1; m < 13; m++)
	{
		sw = 1 + (6.0 * sin((PI * m) / 12.0));
		c[m] = c[m] * sw;
	}
	for (int i = 0; i < 13; i++)
	{
		a_lpc[i] = 0;
		r[i] = 0;
	}
}

void clear()
{
	for (int i = 0; i < 13; i++)
	{
		a_lpc[i] = 0;
		r[i] = 0;
		c[i] = 0;
		for (int j = 0; j < 13; j++)
		{
			alpha[i][j] = 0;
		}
	}
}

void universe(char filename[])
{
	int count = 0;
	FILE *in;
	in = fopen(filename, "r");
	if (!in)
		printf("Unable to access %s\n", filename);

	for (int i = 0; i < MAX; i++)
	{
		fscanf(in, "%Lf", &n_amp[i]);
	}

	normalise_Amplitude();
	hammingWindow();

	for (int l = 0; l < 13; l++)
	{
		for (int p = 0; p < MAX - l; p++)
		{
			r[l] = r[l] + ((n_amp[p]) * n_amp[p + l]);
		}
	}

	durbinAlgo();
	calculateCepstal();
	for (int i = 1; i < 13; i++)
	{
		fprintf(uni, "%Lf,", c[i]);
		universeIn[line][i - 1] = c[i];
	}
	line++;
	fprintf(uni, "\n");
	clear();

	while (!feof(in))
	{
		for (int i = 0; i < 240; i++)
		{
			n_amp[i] = n_amp[80 + i];
		}

		for (int i = 0; i < 80; i++)
		{
			if (!feof(in))
			{
				fscanf(in, "%Lf", &n_amp[i + 240]);
				count++;
			}
		}
		if (count == 80)
		{
			normalise_Amplitude();
			hammingWindow();

			for (int l = 0; l < 13; l++)
			{
				for (int p = 0; p < MAX - l; p++)
				{
					r[l] = r[l] + ((n_amp[p]) * n_amp[p + l]);
				}
			}

			durbinAlgo();
			calculateCepstal();
			for (int i = 1; i < 13; i++)
			{
				fprintf(uni, "%Lf,", c[i]);
				universeIn[line][i - 1] = c[i];
			}
			line++;
			fprintf(uni, "\n");
			clear();
			count = 0;
		}
	}

	fclose(in);
}

void init()
{
	calculate_DCShift("DC_Shift.txt");
	for (int i = 0; i < 10; i++)
	{
		for (int j = 0; j < 10; j++)
		{
			inputFile1[29] = char(48 + i);
			inputFile1[31] = char(48 + j);
			universe(inputFile1);
		}
		for (int j = 0; j < 10; j++)
		{

			inputFile2[29] = char(48 + i);
			inputFile2[32] = char(48 + j);
			universe(inputFile2);
		}
	}
}

//@brief - finds the index of min tokura distance
// return - index of min tokura distance
int minIndex(int k)
{
	int minIndex = 0;
	double min = distance[0];
	for (int i = 1; i < k; i++)
	{
		if (distance[i] < min)
		{
			minIndex = i;
			min = distance[i];
		}
	}
	return minIndex;
}

//@brief - calculate tokura distance with codebook
//@parma - a - vector under consideration
void tokuraDist(long double arr[], int k)
{
	long double t;
	int index;
	for (int i = 0; i < k; i++)
	{
		for (int j = 0; j < 12; j++)
		{
			t = (arr[j] - centroid[i][j]) * (arr[j] - centroid[i][j]);
			distance[i] = distance[i] + tokura_wt[j] * t;
		}
	}
}

//@brief - classifies the vectors in the univers (size n), into k buckets
//@params - k - size of codebook
//@params - n - size of universe
int K_means(int k, int n)
{
	int m = 0, index;
	long double dist_sum = 0;

	for (int i = 0; i < 1000; i++)
	{
		dist_hist[i] = 0;
	}

	for (int i = 0; i < k; i++)
	{
		distortion[i] = 0;
		v_count[i] = 0;
	}

	for (int i = 0; i < n; i++)
	{
		tokuraDist(universeIn[i], k);
		index = minIndex(k);

		for (int j = 0; j < 12; j++)
		{
			codebook[index][j] += universeIn[i][j];
		}

		v_count[index]++;

		distortion[index] += distance[index];

		for (int j = 0; j < k; j++)
		{
			distance[j] = 0;
		}
	}

	for (int i = 0; i < k; i++)
	{
		dist_sum += distortion[i];
	}

	dist_hist[m] = dist_sum / n;
	dist_sum = 0;
	for (int i = 0; i < k; i++)
	{
		for (int j = 0; j < 12; j++)
		{
			centroid[i][j] = codebook[i][j] / v_count[i];
		}
	}

	for (int i = 0; i < k; i++)
	{
		distortion[i] = 0;
		v_count[i] = 0;
	}

	for (int i = 0; i < k; i++)
	{
		for (int j = 0; j < 12; j++)
		{
			codebook[i][j] = 0;
		}
	}

	do
	{
		m++;
		for (int i = 0; i < n; i++)
		{
			tokuraDist(universeIn[i], k);
			index = minIndex(k);

			for (int j = 0; j < 12; j++)
			{
				codebook[index][j] += universeIn[i][j];
			}

			v_count[index]++;

			distortion[index] += distance[index];

			for (int j = 0; j < k; j++)
			{
				distance[j] = 0;
			}
		}
		for (int i = 0; i < k; i++)
		{
			dist_sum += distortion[i];
		}

		dist_hist[m] = dist_sum / n;
		dist_sum = 0;

		for (int i = 0; i < k; i++)
		{
			for (int j = 0; j < 12; j++)
			{
				centroid[i][j] = codebook[i][j] / v_count[i];
			}
		}

		for (int i = 0; i < k; i++)
		{
			distortion[i] = 0;
			v_count[i] = 0;
		}

		for (int i = 0; i < k; i++)
		{
			for (int j = 0; j < 12; j++)
			{
				codebook[i][j] = 0;
			}
		}
	} while (fabs(dist_hist[m] - dist_hist[m - 1]) > delta);

	return 1;
}

//@brief - sets intial codebook
//@params - n - size of universe
void intialisation(int n)
{

	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < 12; j++)
		{
			centroid[l_index][j] += universeIn[i][j];
		}
	}

	for (int i = 0; i < 12; i++)
	{
		centroid[l_index][i] /= n;
	}
}

int LBG(int n, int k)
{
	intialisation(n);

	int x = 1;

	while (x < k)
	{

		// splits codebook
		for (int i = 0; i < x; i++)
		{
			l_index++;
			for (int j = 0; j < 12; j++)
			{
				centroid[l_index][j] = centroid[i][j] * (1 + epsilon);
			}
			for (int j = 0; j < 12; j++)
			{
				centroid[i][j] = centroid[i][j] * (1 - epsilon);
			}
		}
		x *= 2;

		K_means(x, n);
	}
	FILE *cb;
	cb = fopen("Codebook.csv", "w+");
	if (!cb)
		printf("Unable to access Codebook\n");

	for (int i = 0; i < 32; i++)
	{
		for (int j = 0; j < 12; j++)
		{
			fprintf(cb, "%Le, ", centroid[i][j]);
		}
		fprintf(cb, "\n");
	}
	fclose(cb);
	return 1;
}

// find the index of min distance
int min()
{
	int minIndex = 0;
	long double min = distance[0];
	for (int i = 1; i < 32; i++)
	{
		if (distance[i] < min)
		{
			minIndex = i;
			min = distance[i];
		}
	}
	return minIndex;
}

void readAB()
{
	FILE *fp;
	fp = fopen("intialModel/a.txt", "r");

	for (int i = 1; i <= N; i++)
	{
		for (int j = 1; j <= N; j++)
		{
			fscanf(fp, "%Le", &a[i][j]);
		}
	}
	fclose(fp);
	fp = fopen("intialModel/b.txt", "r");
	for (int i = 1; i <= N; i++)
	{
		for (int j = 1; j <= M; j++)
		{
			fscanf(fp, "%Le", &b[i][j]);
		}
	}
	fclose(fp);
	fp = fopen("intialModel/pi.txt", "r");
	for (int i = 1; i <= N; i++)
	{
		fscanf(fp, "%Le", &pi[i]);
	}
	fclose(fp);
}

void findObsSeq(FILE *fp, FILE *temp)
{
	int count = 0;
	for (int i = 0; i < MAX; i++)
	{
		fscanf(fp, "%Lf", &n_amp[i]);
	}

	normalise_Amplitude();
	hammingWindow();
	for (int l = 0; l < 13; l++)
	{
		for (int p = 0; p < MAX - l; p++)
		{
			r[l] = r[l] + ((n_amp[p]) * n_amp[p + l]);
		}
	}

	durbinAlgo();
	calculateCepstal();
	tokuraDist(c, 32);
	o[T++] = min() + 1;
	clear();
	for (int i = 0; i < 32; i++)
	{
		distance[i] = 0;
	}

	while (!feof(fp))
	{
		for (int i = 0; i < 240; i++)
		{
			n_amp[i] = n_amp[80 + i];
		}

		for (int i = 0; i < 80; i++)
		{
			if (!feof(fp))
			{
				fscanf(fp, "%Lf", &n_amp[i + 240]);
				count++;
			}
		}
		if (count == 80)
		{
			normalise_Amplitude();
			hammingWindow();

			for (int l = 0; l < 13; l++)
			{
				for (int p = 0; p < MAX - l; p++)
				{
					r[l] = r[l] + ((n_amp[p]) * n_amp[p + l]);
				}
			}

			durbinAlgo();
			calculateCepstal();
			tokuraDist(c, 32);
			o[T++] = min() + 1;
			clear();
			for (int i = 0; i < 32; i++)
			{
				distance[i] = 0;
			}
			count = 0;
		}
	}

	fprintf(logfp, "Intial observation Sequence\n");
	for (int i = 1; i < T; i++)
	{
		fprintf(logfp, "%d\t", o[i]);
		fprintf(temp, "%d\t", o[i]);
	}
	fprintf(logfp, "\n");
	fprintf(temp, "\n");
}

long double forwardProcedure(FILE *fp)
{
	long double sum = 0;

	for (int i = 1; i <= N; i++)
	{
		alphaf[1][i] = pi[i] * b[i][o[1]];
	}

	for (int t = 1; t <= T - 1; t++)
	{
		for (int j = 1; j <= N; j++)
		{
			for (int i = 1; i <= N; i++)
			{
				sum += (alphaf[t][i] * a[i][j]);
			}
			alphaf[t + 1][j] = sum * b[j][o[t + 1]];
			sum = 0;
		}
	}

	sum = 0;

	for (int m = 1; m <= N; m++)
	{
		sum += alphaf[T][m];
	}
	fprintf(logfp, "\nP(O|lamda) = %Le\n\n", sum);
	fprintf(fp, "\nP(O|lamda) = %Le\n\n", sum);
	for (int i = 1; i <= T; i++)
	{
		for (int j = 1; j <= N; j++)
		{
			fprintf(fp, "%Le ", alphaf[i][j]);
		}
		fprintf(fp, "\n");
	}
	fprintf(fp, "\n\n");

	return sum;
}

void backwardProcedure(FILE *fp)
{
	long double sum = 0;
	for (int i = 1; i <= N; i++)
	{
		beta[T][i] = 1;
	}

	for (int t = T - 1; t >= 1; t--)
	{
		for (int i = 1; i <= N; i++)
		{
			for (int j = 1; j <= N; j++)
			{
				sum += (a[i][j] * b[j][o[t + 1]] * beta[t + 1][j]);
			}
			beta[t][i] = sum;
			sum = 0;
		}
	}

	for (int i = 1; i <= T; i++)
	{
		for (int j = 1; j <= N; j++)
		{
			fprintf(fp, "%Le ", beta[i][j]);
		}
		fprintf(fp, "\n");
	}
	fprintf(fp, "\n\n");
}

void calculateGamma(FILE *fp)
{
	long double temp, sum;
	for (int t = 1; t <= T; t++)
	{
		sum = 0;
		for (int k = 1; k <= N; k++)
		{
			sum += (alphaf[t][k] * beta[t][k]);
		}
		for (int i = 1; i <= N; i++)
		{
			temp = alpha[t][i] * beta[t][i];

			gamma[t][i] = (alphaf[t][i] * beta[t][i]) / sum;
		}
	}

	for (int i = 1; i <= T; i++)
	{
		for (int j = 1; j <= N; j++)
		{
			fprintf(fp, "%Le ", gamma[i][j]);
		}
		fprintf(fp, "\n");
	}
	fprintf(fp, "\n\n");
}

void viterbi(FILE *fp)
{
	long double max, temp;
	int maxIndex;

	for (int i = 1; i <= N; i++)
	{
		del[1][i] = pi[i] * b[i][o[1]];
		si[1][i] = 0;
	}

	for (int t = 2; t <= T; t++)
	{
		for (int j = 1; j <= N; j++)
		{
			max = 0;
			for (int i = 1; i <= N; i++)
			{
				temp = del[t - 1][i] * a[i][j];
				if (max < temp)
				{
					max = temp;
					si[t][j] = i;
				}
			}
			del[t][j] = max * b[j][o[t]];
		}
	}

	max = del[T][1];
	maxIndex = 1;

	for (int i = 2; i <= N; i++)
	{
		if (max < del[T][i])
		{
			max = del[T][i];
			maxIndex = i;
		}
	}

	p_star = max;
	q_star[T] = maxIndex;

	for (int t = T - 1; t >= 1; t--)
	{
		q_star[t] = si[t + 1][q_star[t + 1]];
	}

	fprintf(logfp, "\nQ star");
	fprintf(logfp, "------------\n");
	for (int t = 1; t <= T; t++)
	{
		fprintf(logfp, "%d  ", q_star[t]);
	}
	fprintf(logfp, "\n------------\n");

	fprintf(logfp, "\nP star");
	fprintf(logfp, "\n-----------\n");

	fprintf(logfp, "%Le\n", p_star);

	fprintf(fp, "\nQ star\n");
	fprintf(fp, "\n------------\n");
	for (int t = 1; t <= T; t++)
	{
		fprintf(fp, "%d  ", q_star[t]);
	}
	fprintf(fp, "\n------------\n");

	fprintf(fp, "\nP star\n");
	fprintf(fp, "\n-----------\n");

	fprintf(fp, "%Le\n", p_star);

	p_new = p_star;
}

void clearHMM()
{
	for (int i = 0; i < T_MAX; i++)
	{
		for (int j = 0; j <= N; j++)
		{
			alphaf[i][j] = 0;
			beta[i][j] = 0;
			gamma[i][j] = 0;
			del[i][j] = 0;
			si[i][j] = 0;
			theta[i][j] = 0;
		}
	}
	p_old = p_new = 0;

	for (int i = 0; i < T_MAX; i++)
	{
		for (int j = 0; j <= N; j++)
		{
			for (int k = 0; k < N; k++)
			{
				zeta[i][j][k] = 0;
			}
		}
	}
}

void restimation()
{
	long double sum;

	for (int t = 1; t <= T - 1; t++)
	{
		sum = 0;
		for (int k = 1; k <= N; k++)
		{
			for (int m = 1; m <= N; m++)
			{
				sum += (alphaf[t][k] * a[k][m] * b[m][o[t + 1]] * beta[t + 1][m]);
			}
		}

		for (int i = 1; i <= N; i++)
		{
			for (int j = 1; j <= N; j++)
			{
				if (sum == 0)
					zeta[t][i][j] = 0;
				else
				{
					zeta[t][i][j] = ((alphaf[t][i] * a[i][j] * b[j][o[t + 1]] * beta[t + 1][j]) / sum);
				}
			}
		}
	}

	long double num, den;

	for (int i = 1; i <= N; i++)
	{
		for (int j = 1; j <= N; j++)
		{
			num = den = 0;
			for (int t = 1; t <= T - 1; t++)
			{
				num += zeta[t][i][j];
				den += gamma[t][i];
			}
			if (den == 0)
				a[i][j] = 0;
			else
			{
				a[i][j] = num / den;
			}
		}
	}

	for (int j = 1; j <= N; j++)
	{
		for (int k = 1; k <= M; k++)
		{
			num = den = 0;
			for (int t = 1; t <= T; t++)
			{
				if (q_star[t] == j && o[t] == k)
				{
					num += gamma[t][j];
				}
				den += gamma[t][j];
			}
			if (den == 0)
				b[j][k] = 0;
			else
			{
				b[j][k] = (num / den);
			}
		}
	}
}

void reAdjustment(FILE *fp)
{
	long double max = 0, sum = 0;
	int index = 0;
	for (int i = 1; i <= N; i++)
	{
		for (int j = 1; j <= N; j++)
		{
			sum += a[i][j];
			if (a[i][j] > max)
			{
				max = a[i][j];
				index = j;
			}
		}
		if (sum <= 1)
		{
			a[i][index] += (1 - sum);
		}
		else
		{
			a[i][index] -= (sum - 1);
		}
		sum = max = 0;
		index = 0;
	}

	index = 0;
	max = 0;
	sum = 0;
	for (int i = 1; i <= N; i++)
	{
		for (int j = 1; j <= M; j++)
		{
			if (b[i][j] < MIN)
			{
				b[i][j] = MIN;
			}
			if (max < b[i][j])
			{
				max = b[i][j];
				index = j;
			}
			sum += b[i][j];
		}
		if (sum > 1)
		{
			b[i][index] -= (sum - 1);
		}
		else
		{
			b[i][index] += (1 - sum);
		}
		sum = index = max = 0;
	}
	fprintf(fp, "\nA Martrix\n\n");
	for (int i = 1; i <= N; i++)
	{
		for (int j = 1; j <= N; j++)
		{
			fprintf(fp, "%Le ", a[i][j]);
		}
		fprintf(fp, "\n");
	}

	fprintf(fp, "\nB Martrix\n\n");
	for (int i = 1; i <= N; i++)
	{
		for (int j = 1; j <= M; j++)
		{
			fprintf(fp, "%Le ", b[i][j]);
		}
		fprintf(fp, "\n");
	}
}

void sumA()
{
	for (int i = 1; i <= N; i++)
	{
		for (int j = 1; j <= N; j++)
		{
			a_dash[i][j] += a[i][j];
		}
	}
}

void sumB()
{
	for (int i = 1; i <= N; i++)
	{
		for (int j = 1; j <= M; j++)
		{
			b_dash[i][j] += b[i][j];
		}
	}
}

void stochastic()
{
	long double max = 0, sum = 0;
	int index = 0;
	for (int i = 1; i <= N; i++)
	{
		for (int j = 1; j <= N; j++)
		{
			sum += a_dash[i][j];
			if (a_dash[i][j] > max)
			{
				max = a_dash[i][j];
				index = j;
			}
		}
		if (sum <= 1)
		{
			a_dash[i][index] += (1 - sum);
		}
		else
		{
			a_dash[i][index] -= (sum - 1);
		}
		sum = max = 0;
		index = 0;
	}

	index = 0;
	max = 0;
	sum = 0;
	for (int i = 1; i <= N; i++)
	{
		for (int j = 1; j <= M; j++)
		{
			if (b_dash[i][j] < MIN)
			{
				b_dash[i][j] = MIN;
			}
			if (max < b_dash[i][j])
			{
				max = b_dash[i][j];
				index = j;
			}
			sum += b_dash[i][j];
		}
		if (sum > 1)
		{
			b_dash[i][index] -= (sum - 1);
		}
		else
		{
			b_dash[i][index] += (1 - sum);
		}
		sum = index = max = 0;
	}
}

void training()
{
	FILE *in, *temp, *fop1, *fop2;
	T = 1;
	int iteration = 1;
	int run = 0;
	for (int i = 0; i < 10; i++)
	{
		for (int j = 0; j < 10; j++)
		{
			inputFile1[29] = char(48 + i);
			inputFile1[31] = char(48 + j);
			tempFile1[17] = char(48 + i);
			tempFile1[19] = char(48 + j);

			fprintf(logfp, "Training %s\n", inputFile1);
			fprintf(logfp, "=================\n");

			in = fopen(inputFile1, "r");
			if (!in)
				printf("Unable to access %s\n", inputFile1);

			temp = fopen(tempFile1, "w+");
			if (!temp)
				printf("Unable to access %s\n", tempFile1);

			findObsSeq(in, temp);

			T--;
			readAB();
			fprintf(temp, "Iteraion %d", iteration++);
			fprintf(temp, "\n-----------------\n");
			forwardProcedure(temp);
			backwardProcedure(temp);
			calculateGamma(temp);
			viterbi(temp);
			restimation();
			reAdjustment(temp);
			while (p_old < p_new)
			{
				p_old = p_new;
				fprintf(temp, "Iteraion %d", iteration++);
				fprintf(temp, "\n-----------------\n");
				forwardProcedure(temp);
				backwardProcedure(temp);
				calculateGamma(temp);
				viterbi(temp);
				restimation();
				reAdjustment(temp);
			}
			T = 1;
			iteration = 1;
			sumA();
			sumB();
			clearHMM();
			fprintf(temp, "=================\n");
			fprintf(logfp, "=================\n");
			run = 0;
			fclose(temp);
			fclose(in);
		}
		for (int j = 0; j < 10; j++)
		{
			inputFile2[29] = char(48 + i);
			inputFile2[32] = char(48 + j);
			tempFile2[17] = char(48 + i);
			tempFile2[20] = char(48 + j);
			fprintf(logfp, "Training %s\n", inputFile2);
			fprintf(logfp, "=================\n");

			in = fopen(inputFile2, "r");
			if (!in)
				printf("Unable to access %s\n", inputFile2);

			temp = fopen(tempFile2, "w+");
			if (!temp)
				printf("Unable to access %s\n", tempFile2);

			findObsSeq(in, temp);
			T--;

			readAB();
			fprintf(temp, "Iteraion %d", iteration++);
			fprintf(temp, "\n-----------------\n");
			forwardProcedure(temp);
			backwardProcedure(temp);
			calculateGamma(temp);
			viterbi(temp);
			restimation();
			reAdjustment(temp);
			while (p_old < p_new)
			{
				p_old = p_new;
				fprintf(temp, "Iteraion %d", iteration++);
				fprintf(temp, "\n-----------------\n");
				forwardProcedure(temp);
				backwardProcedure(temp);
				calculateGamma(temp);
				viterbi(temp);
				restimation();
				reAdjustment(temp);
			}
			T = 1;
			iteration = 1;
			fprintf(temp, "=================\n");
			fprintf(logfp, "=================\n");
			sumA();
			sumB();
			run = 0;
			clearHMM();
			fclose(temp);
			fclose(in);
		}
		finalFile1[23] = char(48 + i);
		finalFile2[23] = char(48 + i);
		fop1 = fopen(finalFile1, "w+");
		if (!fop1)
			printf("Unable to access %s\n", finalFile1);
		fop2 = fopen(finalFile2, "w+");
		if (!fop2)
			printf("Unable to access %s\n", finalFile2);

		for (int k = 1; k <= N; k++)
		{
			for (int l = 1; l <= N; l++)
			{
				a_dash[k][l] /= 20;
			}
		}
		for (int k = 1; k <= N; k++)
		{
			for (int l = 1; l <= M; l++)
			{
				b_dash[k][l] /= 20;
			}
		}
		stochastic();
		for (int k = 1; k <= N; k++)
		{
			for (int l = 1; l <= N; l++)
			{
				fprintf(fop1, "%Le ", a_dash[k][l]);
				a_dash[k][l] = 0;
			}
			fprintf(fop1, "\n");
		}
		for (int k = 1; k <= N; k++)
		{
			for (int l = 1; l <= M; l++)
			{
				fprintf(fop2, "%Le ", b_dash[k][l]);
			}
			fprintf(fop2, "\n");
		}
		fclose(fop1);
		fclose(fop2);
	}
}

void readCodebook(FILE *cb)
{
	for (int i = 0; i < 32; i++)
	{
		for (int j = 0; j < 12; j++)
		{
			fscanf(cb, "%Le,", &centroid[i][j]);
		}
	}
}

void readModel(int digit)
{
	FILE *afp, *bfp;
	long double temp;
	finalFile1[23] = char(48 + digit);
	finalFile2[23] = char(48 + digit);

	afp = fopen(finalFile1, "r");
	if (!afp)
		printf("Unable to access %s\n", finalFile1);

	bfp = fopen(finalFile2, "r");
	if (!bfp)
		printf("Unable to access %s\n", finalFile2);

	for (int i = 1; i <= N; i++)
	{
		for (int j = 1; j <= N; j++)
		{
			fscanf(afp, "%Le", &a[i][j]);
		}
	}

	for (int i = 1; i <= N; i++)
	{
		for (int j = 1; j <= M; j++)
		{
			fscanf(bfp, "%Le", &b[i][j]);
		}
	}

	fclose(afp);
	fclose(bfp);
}

int maxPro(FILE *fp)
{
	int index = 0;
	long double maxValue = pTest[0];
	pTest[0] = 0.0;
	for (int i = 1; i < 10; i++)
	{
		if (maxValue < pTest[i])
		{
			index = i;
			maxValue = pTest[i];
		}
		pTest[i] = 0;
	}
	fprintf(fp, "\nMax probability = %Le\n", maxValue);
	return index;
}

void clrTest()
{
	for (int i = 1; i <= N; i++)
	{
		for (int j = 1; j <= N; j++)
		{
			a[i][j] = 0;
			;
		}
	}

	for (int i = 1; i <= N; i++)
	{
		for (int j = 1; j <= M; j++)
		{
			b[i][j] = 0;
		}
	}

	for (int i = 0; i < T_MAX; i++)
	{
		o[i] = 0;
	}
}

void testing()
{
	FILE *testfp, *temp;
	int max;
	int accuracy = 0;
	for (int i = 0; i < 10; i++)
	{
		for (int j = 0; j < 10; j++)
		{
			T = 1;
			test[29] = char(48 + i);
			test[32] = char(48 + j);
			testfp = fopen(test, "r");
			if (!testfp)
				printf("Unable to access %s\n", test);

			tempFile3[17] = char(48 + i);
			tempFile3[20] = char(48 + j);
			temp = fopen(tempFile3, "w+");
			if (!temp)
				printf("Unable to access %s\n", tempFile3);

			findObsSeq(testfp, temp);
			T--;
			for (int k = 0; k < 10; k++)
			{
				readModel(k);
				pTest[k] = forwardProcedure(temp);
			}
			max = maxPro(temp);
			fprintf(temp, "\nUtterance recognised as %d", max);
			if (max == i)
				accuracy++;
			clrTest();
			T = 1;
			fclose(temp);
			fclose(testfp);
		}
	}
	printf("Accuracy is %d %\n", accuracy);
}

void liveTesting()
{
	FILE *finput, *output, *rfp, *temp;
	int amp;
	int max;
	int count = 0;

	system("Recording_Module.exe 3 input_file.wav input_file.txt");

	finput = fopen("input_file.txt", "r+");
	output = fopen("recording_output.txt", "w+");

	while (!feof(finput))
	{
		fscanf(finput, "%d", &amp);
		if (amp > 200)
			break;
	}

	while (!feof(finput))
	{
		fscanf(finput, "%d", &amp);
		if (count > 0)
		{
			if (amp < 100)
			{
				count++;
			}
			else
			{
				count = 0;
			}
		}
		else
		{
			if (amp < 100)
			{
				count++;
			}
			else
			{
				fprintf(output, "%d\n", amp);
			}
		}
		if (count >= 100)
		{
			break;
		}
	}

	fclose(finput);
	fclose(output);

	rfp = fopen("recording_output.txt", "r+");
	temp = fopen("temp/Live_Record_Log.txt", "w+");
	T = 1;
	findObsSeq(rfp, temp);
	T--;
	for (int k = 0; k < 10; k++)
	{
		readModel(k);
		pTest[k] = forwardProcedure(temp);
	}
	max = maxPro(temp);
	fprintf(temp, "\nUtterance recognised as %d", max);
	clrTest();
	T = 1;
	fclose(temp);
	fclose(rfp);
}

int _tmain(int argc, _TCHAR *argv[])
{
	FILE *cbin;
	int choice;

	// // Training
	// cbin = fopen("Codebook.csv", "r");
	// if (!cbin)
	// 	printf("Unable to access Codebook\n");

	// readCodebook(cbin);

	// training();

	// fclose(cbin);

	printf("DIGIT RECOGNITION SYSTEM\n");
	printf("==========================\n");
	printf("\nChoose a option\n");
	printf("\n1. Create Codebook");
	printf("\n2. Train Model");
	printf("\n3. Test Pre-recorded file");
	printf("\n4. Test with live recording");
	printf("\n5. Exit\n\n");
	scanf("%d", &choice);

	logfp = fopen("logfp.txt", "w+");
	if (!logfp)
		printf("Unable to access logfp\n");
	printf("\nLog File Created\n");

	while (1)
	{
		switch (choice)
		{
		case 1:
			uni = fopen("Universe.csv", "w+");
			if (!uni)
				printf("Unable to access Universe\n");

			init();
			printf("Universe Created\n");
			fclose(uni);
			LBG(9106, 32);
			printf("Codebook Created\n");
			exit(0);
			break;
		case 2:
			cbin = fopen("Codebook.csv", "r");
			if (!cbin)
				printf("Unable to access Codebook\n");
			else
				printf("Reading from codebook\n");
			readCodebook(cbin);
			printf("Training Started\n");
			training();
			printf("Training Ended\nModel stored in Final Model Folder");
			fclose(cbin);
			exit(0);
			break;
		case 3:
			cbin = fopen("Codebook.csv", "r");
			if (!cbin)
				printf("Unable to access Codebook\n");
			else
				printf("Reading from codebook\n");
			readCodebook(cbin);
			printf("Training Started\n");
			training();
			printf("Training Ended\nModel stored in Final Model Folder\n");
			testing();
			printf("Testing Ended\nOutput stored in Temp Folder\n");
			getch();
			fclose(cbin);
			exit(0);
			break;
		case 4:
			cbin = fopen("Codebook.csv", "r");
			if (!cbin)
				printf("Unable to access Codebook\n");
			else
				printf("Reading from codebook\n");
			readCodebook(cbin);
			printf("Training Started\n");
			training();
			printf("Training Ended\nModel stored in Final Model Folder\n");
			liveTesting();
			fclose(cbin);
			exit(0);
			break;
		case 5:
			fclose(logfp);
			exit(0);
			break;
		default:
			printf("Enter Valid Number\n");
			break;
		}
	}
	printf("End of Program\n");
	fclose(logfp);
	getch();
	return 0;
}
