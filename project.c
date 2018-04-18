// MUSI4633 - Final Project - Phase Distortion Synthesis
// Professor: Georg Boenn
// By: Alex Hochheiden & Keegan Herperger
// Spring 2018

#include <math.h>
#include <stdio.h> 
#include <stdlib.h>
#include <string.h>
#include "wave.h"  /* header file containing the wavehead struc 
					  and the declaration of update_header() */
FILE* fpout;
const double TWO_PI = 2 * acos(-1);

void envelope(double *arr, int sr, int dur);
//void carriergen(int end, int blockframes, short* audioblock, int freq, int sr, int dur, double* env, int bd);
//void modgen(int end, int blockframes, short* audioblock, int freq, int sr, int dur, double* env, int bd);
//void breaks(float* pointList, float* levelList, int length);
//void addsyn(int dur);

int main(int argc, char** argv)
{
	short  *audioblock;      /* audio memory pointer */
	int     end, i, j, n;       /* dur in frames, counter vars */
	int     sr = 44100;      /* sampling rate */
	int     blockframes = 256; /* audio block size in frames */
	int     databytes;  /* audio data in bytes */
	unsigned int ndx = 0, adsrIndex = 0;   /* phase index for synthesis */
	float   dur, freq, phasorFreq, sustain, prevBaseFrequencyCounter; /* duration, frequency, phasor frequency, sustain */
	int attack, decay, release, depth, shape;   /* attack, decay, release, bit-depth, shape*/

	wavehead *header;

	if(argc != 6)
	{
		printf("usage: %s outfile dur freq phasorFreq bit-depth samplerate\n", argv[0]);
		exit(-1);
	}

	fpout = fopen(argv[1], "wb");
	dur = atof(argv[2]);
	freq = atof(argv[3]);
	phasorFreq = atof(argv[4]);
	depth = atoi(argv[5]);
	sr = atof(argv[6]);
	end = (int)(dur*sr);
	audioblock = (short *)malloc(sizeof(short)*blockframes);
	prevBaseFrequencyCounter = -999;

	/* set the data size */
	databytes = end * sizeof(short);

	//allocate space for the envelope
	double* env = (double*)malloc(end * sizeof(double));

	/* write the header */
	header = (wavehead*)malloc(sizeof(wavehead));
	update_header(header, sr, 1, 16, databytes);
	fwrite(header, 1, sizeof(wavehead), fpout);

	// adsr
	envelope(env, sr, dur);

	for(i = 0; i < end; i += blockframes)
	{
		for(j = 0; j < blockframes; j++, ndx++, adsrIndex++)
		{
			float baseFrequencyCounter = 0;
			audioblock[j] = 0;
			
			for(int n = 1; n < 100; n ++)
			{
				if(n * freq >= sr / 2)	
				{ 
					// nyquist frequency check
					break;
				}

				// both sawtooth?
				baseFrequencyCounter += 16000 * pow(-1, (n + 1)) / n * sin(ndx * TWO_PI * n * freq / sr);
				audioblock[j] = 16000 * pow(-1, (n + 1)) / n * sin(ndx * TWO_PI * n * phasorFreq / sr);
			}

			if(ndx >= dur*sr)
			{
				break;
			}
			else
			{
				audioblock[j] *= (-baseFrequencyCounter) * env[adsrIndex];

				// make sure this works, should reset after the drop at the end of a sawtooth wave
				if(baseFrequencyCounter > prevBaseFrequencyCounter)
				{
					prevBaseFrequencyCounter = baseFrequencyCounter;
					continue;
				}
				else
				{
					// reset phase on end of base frequency period
					ndx = 0;
					prevBaseFrequencyCounter = -999;
				}				 
			}
		}
		fwrite(audioblock, sizeof(short), blockframes, fpout);
	}

	free(audioblock);
	free(header);
	fclose(fpout);

	return 0;
}

void envelope(double *arr, int sr, int dur)
{
	int a, d, r;
	double s;

	printf("Attack time in ms:");
	scanf("%d", &a);
	printf("Decay time in ms:");
	scanf("%d", &d);
	printf("Sustain level (0-1):");
	scanf("%lf", &s);
	printf("Release time in ms:");
	scanf("%d", &r);

	int period = sr / 1000;
	if(a + d + r >= dur * 1000)
	{
		printf("Invalid envelope length. Amplitude will be constant. \n");
		int i = 0;
		while(arr[i])
		{
			arr[i] = 0.8;
			i++;
		}
	}
	else
	{
		double inc = 1. / (period*a);
		for(int i = 0; i < period*a; i++)
		{
			arr[i] = inc*i;
		}

		inc = (1. - s) / (period*d);
		for(int i = period*(a + d); i > period*a; i--)
		{
			arr[i] = 1. - (inc*i);
		}
		for(int i = period*d; i < (sr*dur - period*r); i++)
		{
			arr[i] = s;
		}

		inc = s / (period*r);
		for(int i = period*r; i > 0; i--)
		{
			arr[sr*dur - i] = inc*i;
		}
	}
}

/*
void modgen(int end, int blockframes, short* audioblock, int freq, int sr, int dur, double* env, int bd)
{
	int jnum;
	printf("How many jumps would you like in the modulator signal?");
	scanf("%i", &jnum);

	int spc = sr*freq; //samples per cycle
	int cycleLength = sr*dur;

	float* wave = (float*)malloc(spc * sizeof(float));		//array of one phasor cycle
	float* pointList = (float*)malloc(jnum * sizeof(float));
	float* levelList = (float*)malloc(jnum * sizeof(float));

	breaks(pointList, levelList, jnum);
	linsegs(&pointList, &levelList, wave, slope(levelList), cycleLength, spc);

	envelope(env, sr, dur);	// **** We can either keep this as an envelope
							// or replace it with an LFO when we get to that point ***
	// experimenting with new write function here
	for(i = 0; i < end; i += blockframes)
	{
		for(j = 0; j < blockframes; j++, ndx++)
		{
			audioblock[j] = 0;//8000*sin(ndx*TWO_PI*freq/sr);

			for(int n = 1; n < 200; n += 2)
			{ //this is adding sinudoidal harmonics
				if(n*freq >= sr / 2)
				{ //nyquist frequency check
					break;
				}

				audioblock[j] += 16000 * (1. / n)*sin((ndx*TWO_PI*n*freq) / sr);
			}

			if(ndx >= dur*sr)
			{ 
				//this is ending the loop before it goes out of bounds
				continue;
			}
			else
			{
				audioblock[j] *= env[ndx];
			}
		}
		fwrite(audioblock, sizeof(short), blockframes, fpout);
	}
}

void breaks(float* pointList, float* levelList, int length)
//This is the menu for the number of jumps the user wants in the modulator
{
	bool flag = 0;
	int freq;

	printf("Input phasor frequency:");
	scanf("%i", &freq);

	for(int i = 0; i < length; i++){	
		printf("Input next jump point as a percentage of phase (0.-1.):");
		scanf("%f", &pointList[i]);

		printf("Input jump destination as a percentage of amplitude (0.-1.):");
		scanf("%f", &levelList[i]);
	}
}

void addsyn(int dur)
{
	int harmnum, relamp;

	printf("How many harmonics would you like in the carrier?");
	scanf("%d", &harmnum);

	printf("Choose the relative harmonic strength:\n 1. 1/n\n 2. 1/n^2\n 3.sqrt(n) (normalised)\n");
	scanf("%d", &relamp);

	// #TODO fix me
	// carriergen(dur, harmnum, relamp);
}

void carriergen(int end, int blockframes, short* audioblock, int freq, int sr, int dur, double* env, int bd)
{
	int i, j, n, ndx = 0;
	envelope(env, sr, dur);

	for(i = 0; i < end; i += blockframes)
	{
		for(j = 0; j < blockframes; j++, ndx++)
		{
			audioblock[j] = 0;

			for(int h = 1; h < 200; h++)
			{
				n = (2 * h) - 1;
				if(n*freq >= sr / 2)
				{
					break;
				}

				audioblock[j] += 16000 * (1. / n)*sin((ndx*TWO_PI*n*freq) / sr);
			}
			if(ndx >= dur*sr)
			{
				continue;
			}
			else
			{
				audioblock[j] *= env[ndx];
			}
		}
		fwrite(audioblock, sizeof(short), blockframes, fpout);
	}
}
*/

void update_header(wavehead* header, int sr, int channels, int precision, int databytes)
{
	header->magic = (*(long *)RIFF_ID);
	header->len0 = databytes + sizeof(wavehead) - 8;
	header->magic1 = (*(long *)WAVE_ID);
	header->magic2 = (*(long *)FMT_ID);
	header->len = 16;
	header->format = 1;
	header->nchns = (short)channels;
	header->rate = sr;
	header->aver = sr*channels*precision / 8;
	header->nBlockAlign = (short)(channels*precision / 8);
	header->bits = (double)precision;
	header->magic3 = (*(long *)DATA_ID);
	header->datasize = databytes;
}