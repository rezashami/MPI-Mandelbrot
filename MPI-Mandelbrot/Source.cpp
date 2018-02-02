//Sequential MandelBrot Set
//Adapted from source code of Bibek Subedi, http://www.programming-techniques.com/2012/03/computer-graphics-fractals-generation.html
//Modified by Ebrahim Ansari, 01/09/2015
//BMP version is written by Ebrahim Ansari
//MPI version is written by Reza Shami Tanha

#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "mpi.h"
#include <iostream>

typedef struct {
	float x, y;
} Complex;

Complex complexSquare(Complex c) {
	Complex cSq;
	cSq.x = c.x * c.x - c.y * c.y;
	cSq.y = 2 * c.x * c.y;
	return cSq;
}

int iterate(Complex zInit, int maxIter) {
	Complex z = zInit;
	int cnt = 0;
	while ((z.x * z.x + z.y * z.y <= 4) && (cnt < maxIter)) {
		z = complexSquare(z);
		z.x += zInit.x;
		z.y += zInit.y;
		cnt++;
	}
	return cnt;
}


void madelbrot(int nx, int ny, int maxIter, float realMin, float realMax, float imagMin, float imagMax, unsigned char* img) {

	static unsigned char color[3];
	float realInc = (realMax - realMin) / nx;
	float imagInc = (imagMax - imagMin) / ny;

	Complex z;
	int x, y, w, h;
	int cnt;
	for (y = 0, z.y = imagMin; y < ny; y++, z.y += imagInc)
	{
		for (x = 0, z.x = realMin; x < nx; x++, z.x += realInc) {
			cnt = iterate(z, maxIter);
			//using cnt to calculate pixel's colour
			w = x;
			h = (ny - 1) - y;
			if (cnt == maxIter)
			{
				color[0] = 0;
				color[1] = 0;
				color[2] = 0;
			}
			else {
				double c = 3 * log((double)cnt) / log(maxIter - 1.0);
				if (c < 1)
				{
					color[0] = 255 * c;
					color[1] = 0;
					color[2] = 0;
				}
				else if (c < 2)
				{
					color[0] = 255;
					color[1] = 255 * (c - 1);
					color[2] = 0;
				}
				else
				{
					color[0] = 255;
					color[1] = 255;
					color[2] = 255 * (c - 2);
				}
			}
			//putting colours in corresponding pixel (w,h)
			img[(w + h*nx) * 3 + 2] = color[0];
			img[(w + h*nx) * 3 + 1] = color[1];
			img[(w + h*nx) * 3 + 0] = color[2];

		}
	}
}

void embarrassedly_master(int world_size, const int iXmax, const int iYmax, const int iterations, const float realMin, const float realMax)
{
	/*Image pointer to create bmp file*/
	unsigned char *img = NULL;

	/*Pointer to file*/
	FILE * fp;

	/*bmp file nmae*/
	char *filename = "new1.bmp";    //bmp version

									/*Set start time to calculate running time*/
	double t1 = MPI_Wtime();
	/*Divide the imaginary section to number of proccesses.*/
	double incerment_imag = 2/world_size;
	const double myImagMin = -1.00;
	const double myImagMax = -0.50;

	/*File section*/
	int filesize = 54 + 3 * iXmax * iYmax;
	if (img)
		delete img;
	img = new unsigned char[3 * iXmax * iYmax];
	memset(img, 0, sizeof(img));

	/*Divide the img array with this pionner*/
	int array_index = 3240000;
	incerment_imag = -0.5;
	for (int i = 1; i < world_size; i++)
	{
		/*Send the incerment_imag number to i proccess as start imaginary number*/
		MPI_Send(&incerment_imag, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);

		/*Send array to i proccess*/
		MPI_Send(img + array_index, 1620000, MPI_UNSIGNED_CHAR, i, 0, MPI_COMM_WORLD);

		/*Decrease the pionner*/
		array_index -= 1620000;

		/*Increase the incerment_imag*/
		incerment_imag += 0.5;
	}
	/*Execute mandelbrot with proccess infromation*/
	madelbrot(iXmax, 300, iterations, realMin, realMax, myImagMin, myImagMax, img + 4860000);
	/*Set the pionner to defualt value*/
	array_index = 3240000;
	for (int i = 1; i < world_size; i++) {
		MPI_Status st;

		/*Recieve the imag array information*/
		MPI_Recv(img + (array_index), 1620000, MPI_UNSIGNED_CHAR, i, 0, MPI_COMM_WORLD, &st);

		/*Increase the array_index*/
		array_index -= 1620000;
	}
	/*Write section*/
	unsigned char bmpfileheader[14] = { 'B', 'M', 0, 0, 0, 0, 0, 0, 0, 0, 54, 0, 0, 0 };
	unsigned char bmpinfoheader[40] = { 40, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 24, 0 };
	unsigned char bmppad[3] = { 0, 0, 0 };

	bmpfileheader[2] = (unsigned char)(filesize);
	bmpfileheader[3] = (unsigned char)(filesize >> 8);
	bmpfileheader[4] = (unsigned char)(filesize >> 16);
	bmpfileheader[5] = (unsigned char)(filesize >> 24);

	bmpinfoheader[4] = (unsigned char)(iXmax);
	bmpinfoheader[5] = (unsigned char)(iXmax >> 8);
	bmpinfoheader[6] = (unsigned char)(iXmax >> 16);
	bmpinfoheader[7] = (unsigned char)(iXmax >> 24);
	bmpinfoheader[8] = (unsigned char)(iYmax);
	bmpinfoheader[9] = (unsigned char)(iYmax >> 8);
	bmpinfoheader[10] = (unsigned char)(iYmax >> 16);
	bmpinfoheader[11] = (unsigned char)(iYmax >> 24);

	fp = fopen(filename, "wb");
	fwrite(bmpfileheader, 1, 14, fp);
	fwrite(bmpinfoheader, 1, 40, fp);
	//writing pixels in BMP file
	for (int i = 0; i < iYmax; i++)
	{
		fwrite(img + (iXmax *(iYmax - i - 1) * 3), 3, iXmax, fp);
		fwrite(bmppad, 1, (4 - (iXmax * 3) % 4) % 4, fp);
	}
	fclose(fp);

	/*Print the runnig time duretion*/
	std::cout << "The time is:" << MPI_Wtime() - t1 << " from proccess: 0\n";
}

void embarrassedly_slave(int world_size, const int iXmax, const int iYmax, const int iterations, const float realMin, const float realMax)
{
	/*Recieved array stored in this variable*/
	unsigned char* myBuffer = new unsigned char[1620000];

	/*Recieved imaginary information stored in this variables*/
	double myImagMin;
	double myImagMax;

	MPI_Status st2;
	MPI_Status st3;

	/*Recieve the imaginary information */
	MPI_Recv(&myImagMin, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &st2);

	/*Recieve the array information */
	MPI_Recv(myBuffer, 1620000, MPI_UNSIGNED_CHAR, 0, 0, MPI_COMM_WORLD, &st3);

	/*Calculate end of imaginary number*/
	myImagMax = myImagMin + 0.49;

	/*Execute mandelbrot with proccess infromation*/
	madelbrot(iXmax, 300, iterations, realMin, realMax, myImagMin, myImagMax, myBuffer);

	/*Send the calculated array*/
	MPI_Send(myBuffer, 1620000, MPI_UNSIGNED_CHAR, 0, 0, MPI_COMM_WORLD);
}

void oneline_master(int world_size, const int iXmax, const int iYmax, const int iterations, const float realMin, const float realMax) {
	/*Image pointer to create bmp file*/
	unsigned char *img = NULL;

	/*Pointer to file*/
	FILE * fp;

	/*bmp file nmae*/
	char *filename = "new1.bmp";    //bmp version

									/*Set start time to calculate running time*/
	double t1 = MPI_Wtime();
	/*Divide the imaginary section to number of proccesses.*/
	float myImagMin = -1.00;
	float myImagMax;

	/*File section*/
	int filesize = 54 + 3 * iXmax * iYmax;
	if (img)
		delete img;
	img = new unsigned char[3 * iXmax * iYmax];
	memset(img, 0, sizeof(img));

	int array_index = 1619999 - 2700;
	float yIncrement = 0;
	printf("%f \n",yIncrement);
	for (int section = 0; section < 150; section++)
	{
		float temp_Inc = (-1.00 + yIncrement);
		int forward_array = array_index;
		for (int i = 1; i < world_size; i++)
		{
			/*Send the incerment_imag number to i proccess as start imaginary number*/
			MPI_Send(&temp_Inc, 1, MPI_FLOAT, i, 0, MPI_COMM_WORLD);
			printf("%f is sent to %i\n", temp_Inc, i);
			/*Send array to i proccess*/
			MPI_Send(img + (array_index), 2700, MPI_UNSIGNED_CHAR, i, 0, MPI_COMM_WORLD);

			/*Decrease the pionner*/
			array_index -= 2700;

			/*Increase the incerment_imag*/
			yIncrement += (float)2.00 / 600;
			temp_Inc = (-1.00 + yIncrement);
		}

		yIncrement += (float)2.00 / 600;
		temp_Inc = (-1.00 + yIncrement);
		myImagMin = temp_Inc;
		myImagMax = myImagMin + temp_Inc;
		/*Execute mandelbrot with proccess infromation*/
		madelbrot(iXmax, 1, iterations, realMin, realMax, myImagMin, myImagMax, img + array_index);
		/*Set the pionner to defualt value*/
		array_index -= 2700;
		for (int i = 1; i < world_size; i++) {
			MPI_Status st;
			/*Recieve the imag array information*/
			MPI_Recv(img + (forward_array), 2700, MPI_UNSIGNED_CHAR, i, 0, MPI_COMM_WORLD, &st);

			/*Increase the array_index*/
			forward_array -= 2700;
		}
	}
	printf("%f \n", yIncrement);
	printf("%f \n", array_index);
	/*Divide the img array with this pionner*/
	
	/*Write section*/
	unsigned char bmpfileheader[14] = { 'B', 'M', 0, 0, 0, 0, 0, 0, 0, 0, 54, 0, 0, 0 };
	unsigned char bmpinfoheader[40] = { 40, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 24, 0 };
	unsigned char bmppad[3] = { 0, 0, 0 };

	bmpfileheader[2] = (unsigned char)(filesize);
	bmpfileheader[3] = (unsigned char)(filesize >> 8);
	bmpfileheader[4] = (unsigned char)(filesize >> 16);
	bmpfileheader[5] = (unsigned char)(filesize >> 24);

	bmpinfoheader[4] = (unsigned char)(iXmax);
	bmpinfoheader[5] = (unsigned char)(iXmax >> 8);
	bmpinfoheader[6] = (unsigned char)(iXmax >> 16);
	bmpinfoheader[7] = (unsigned char)(iXmax >> 24);
	bmpinfoheader[8] = (unsigned char)(iYmax);
	bmpinfoheader[9] = (unsigned char)(iYmax >> 8);
	bmpinfoheader[10] = (unsigned char)(iYmax >> 16);
	bmpinfoheader[11] = (unsigned char)(iYmax >> 24);

	fp = fopen(filename, "wb");
	fwrite(bmpfileheader, 1, 14, fp);
	fwrite(bmpinfoheader, 1, 40, fp);
	//writing pixels in BMP file
	for (int i = 0; i < iYmax; i++)
	{
		fwrite(img + (iXmax *(iYmax - i - 1) * 3), 3, iXmax, fp);
		fwrite(bmppad, 1, (4 - (iXmax * 3) % 4) % 4, fp);
	}
	fclose(fp);

	/*Print the runnig time duretion*/
	std::cout << "The time is:" << MPI_Wtime() - t1 << " from proccess: 0\n";
}

void oneline_slave(int world_size, const int iXmax, const int iYmax, const int iterations, const float realMin, const float realMax) {
	for (int  i = 0; i < 150; i++)
	{
		/*Recieved array stored in this variable*/
		unsigned char* myBuffer = new unsigned char[2700];

		/*Recieved imaginary information stored in this variables*/
		float myImagMin;
		float myImagMax;

		MPI_Status st2;
		MPI_Status st3;

		/*Recieve the imaginary information */
		MPI_Recv(&myImagMin, 1, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, &st2);

		/*Recieve the array information */
		MPI_Recv(myBuffer, 2700, MPI_UNSIGNED_CHAR, 0, 0, MPI_COMM_WORLD, &st3);

		/*Calculate end of imaginary number*/
		myImagMax = myImagMin + (float)2.00 / 600;

		/*Execute mandelbrot with proccess infromation*/
		madelbrot(iXmax, 1, iterations, realMin, realMax, myImagMin, myImagMax, myBuffer);

		/*Send the calculated array*/
		MPI_Send(myBuffer, 2700, MPI_UNSIGNED_CHAR, 0, 0, MPI_COMM_WORLD);
	}
	
}

int main(int argc, char** argv) {
	//initial setting
	/*Use for keep information about proccesses*/
	int myid;
	int world_size;

	/*Initialize mpi*/
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);

	/*Max iterations number*/
	const int iterations = 1000;

	/*Real number minimum and maximum*/
	const float realMin = -2.0;
	const float realMax = 1.0;

	/*Set the information of bmp file*/
	const int iXmax = 900;
	const int iYmax = 600;

	/*Master section*/
	if (myid == 0)
	{
		//embarrassedly_master(world_size, iXmax, iYmax, iterations, realMin, realMax);
		oneline_master(world_size, iXmax, iYmax, iterations, realMin, realMax);
	}
	else {
		//embarrassedly_slave(world_size, iXmax, iYmax, iterations, realMin, realMax);
		oneline_slave(world_size, iXmax, iYmax, iterations, realMin, realMax);
	}
	/*Finilized the MPI*/
	MPI_Finalize();
	return 0;
}