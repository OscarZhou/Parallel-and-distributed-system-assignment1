// CalculatePI.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"
#include "mpi.h"
#include "math.h"
#include "stdio.h"
#include "stdlib.h"


long calExpo(long x, long n, long p)
{
	long y = 1;
	for (int i = 0; i < n; ++i)
	{
		y = (y*x) % p;
	}
	return y;
}

long calA(long a, int k, long m)
{
	return calExpo(a, k, m);
}

long calC(long a, long c, int k, long m)
{
	long sumf = 0;
	for (int i = k; i >= 0; --i)
	{
		sumf += calExpo(a, i, m);
	}
	long res = (c*(sumf%m)) % m;
	return res;
}

double generateRandomNumWithLeapFrog1(long seed, int pno, int n, int p)
{
	long a = 1664525;
	long m = long(pow(2.0, 32));
	long c = 1013904223;

	double sumf = 0;
	long i_prev = seed;
	long A = calA(a, p, m);
	long C = calC(a, c, p, m);
	for (int i = pno; i < n; i = i + p)
	{
		long i_next = (A*i_prev + C) % m;
		double x_rand = float(i_next) / (m - 1);
		i_prev = i_next;
		double fx = sqrt(1 - x_rand*x_rand);
		sumf += fx;
	}
	double area = sumf / n;
	double pi_calc = 4 * area;
	return pi_calc;
}

double generateRandomNumWithNaiveMethod(int seed, int n)
{
	long a = 1664525;
	long m = long(pow(2.0, 32));
	long c = 1013904223;

	double sumf = 0;
	long i_prev = seed;
	for (int i = 0; i < n;++i)
	{
		long i_next = (a*i_prev + c) % m;
		double x_rand = float(i_next) / (m - 1);
		i_prev = i_next;
		double fx = sqrt(1 - x_rand*x_rand);
		sumf += fx;
	}
	double area = sumf / n;
	double pi_calc = 4 * area;
	return pi_calc;
}


long* generateRandomNumWithLeapFrog(long seed, int n, int p, long arraySeed[100]);

int main(int argc, char* argv[])
{
	double result=0, ret0, ret1;
	MPI_Status Stat;// status variable, so operations can be checked

	MPI_Init(&argc, &argv);

	int world_size;
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);

	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);


	// Firstly, generate the seed for other processors with leapfrog
	fprintf(stdout, "if(rank=0\n");
	if (rank==0)
	{
		// Generate seeds of corresponding numbers with the numbers of processors 
		long arraySeed[100] = {};
		long* basicSeed = generateRandomNumWithLeapFrog(123, world_size, 1, arraySeed);

		// Master sends seed to all the slave processes
		for (int i = 1; i < world_size; i++)
		{
			fprintf(stdout, "seed[%d]=%ld\n", i, basicSeed+i);
			MPI_Send(basicSeed+i, sizeof(long), MPI_LONG, i, 0, MPI_COMM_WORLD);
		}

		// Master calculates its own partial pi
		ret0 = generateRandomNumWithLeapFrog1(basicSeed[0], 0, 100, world_size);

		result += ret0;
		// Receive nodes from all nodes 
		for (int i = 1; i < world_size; i++)
		{
			MPI_Recv(&ret1, sizeof(double), MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &Stat);
			result += ret1;
			fprintf(stdout, "node=%ld: result=%.4f: ret1=%.4f\n", i, result, ret1);

		}
		fprintf(stdout, "The pi calculated by %ld random numbers is %.4f \n", 100, result);
	}
	else // this is not the master
	{
		
		ret1 = 0;
		for (int i = rank; i < world_size;i++)
		{
			long slaveseed;
			MPI_Recv(&slaveseed, sizeof(long), MPI_LONG, 0, 0, MPI_COMM_WORLD, &Stat);
			ret1 = ret1 + generateRandomNumWithLeapFrog1(slaveseed, rank, 100, world_size);

		}
		MPI_Send(&ret1, sizeof(double), MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
	}

	MPI_Finalize();

}



long* generateRandomNumWithLeapFrog(long seed, int n, int p, long arraySeed[100])
{
	long a = 1664525;
	long m = long(pow(2.0, 32));
	long c = 1013904223;

	long sumf = 0;
	long i_prev = seed;
	long A = calA(a, p, m);
	long C = calC(a, c, p, m);

	//long arraySeed[100];
	for (int i = 0; i < n; i = i + p)
	{
		long i_next = (A*i_prev + C) % m;
		arraySeed[i] = i_next;
		/*double x_rand = float(i_next) / (m - 1);
		i_prev = i_next;
		long fx = sqrt(1 - x_rand*x_rand);
		sumf += fx;*/
	}
	//double area = sumf / n;
	//double pi_calc = 4 * area;
	//return pi_calc;
	return arraySeed;
}