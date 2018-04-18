// Course: MUSI4633 
// Project: Final Project - Phase Distortion Synthesis
// Professor: Georg Boenn
// By: Alex Hochheiden & Keegan Herperger
// Semester: Spring 2018

#include <math.h> 
#include <stdio.h> 
#include <stdlib.h>
#include <string.h>
#include "wave.h"  /* header file containing the wavehead struc 
					  and the declaration of update_header() */
FILE* fpout;

void envelope(double *arr, int sr, int dur);
void modulator(double* invertedBaseFreqTable, int sr, int freq, int samplesInBaseFrequencyPeriod, double TWO_PI);
void carrier(double *arr, int sr, int phasorFreq, int samplesInPhasorFrequencyPeriod, double TWO_PI);

int main(int argc, char** argv)
{
	short   *audioblock;      /* audio memory pointer */
	int     end, i, j, n;       /* dur in frames, counter vars */
	int     sr = 44100, samplesInBaseFrequencyPeriod = 0, samplesInPhasorFrequencyPeriod = 0;      /* sampling rate, number of samples in a full cycle at the base frequency */
	int     blockframes = 256; /* audio block size in frames */
	int     databytes;  /* audio data in bytes */
	int     channels = 1;
	unsigned int phaseIndex = 0, adsrIndex = 0;   /* phase index for synthesis, cumulativep hase index used by adsr envelope */
	const double TWO_PI = 2 * acos(-1);
	float   dur, freq, phasorFreq, sustain; /* duration, frequency, phasor frequency, sustain */
	int attack, decay, release, depth = 16;   /* attack, decay, release, bit-depth*/

	wavehead *header;

	if(argc != 6)
	{
		printf("usage: %s outfile dur freq phasorFreq sampleRate\n", argv[0]);
		exit(-1);
	}

	fpout = fopen(argv[1], "wb");
	dur = atof(argv[2]);
	freq = atof(argv[3]);
	phasorFreq = atof(argv[4]);
	sr = atof(argv[5]);
	end = (int)(dur*sr);
	samplesInBaseFrequencyPeriod = (int)(sr/freq);
	samplesInPhasorFrequencyPeriod = (int)(sr/phasorFreq);
	audioblock = (short *)malloc(sizeof(short)*blockframes);

	/* set the data size */
	databytes = end * channels * sizeof(short);

	//allocate space for the envelope
	double* env = (double*)malloc(end * sizeof(double));

	/* write the header */
	header = (wavehead*)malloc(sizeof(wavehead));
	update_header(header, sr, channels, 16, databytes);
	fwrite(header, 1, sizeof(wavehead), fpout);

	// adsr
	envelope(env, sr, dur);

	double* invertedBaseFreqTable = (double*)malloc(samplesInBaseFrequencyPeriod * sizeof(double));
	double* phasorFreqTable = (double*)malloc(samplesInPhasorFrequencyPeriod * sizeof(double));

	modulator(invertedBaseFreqTable, sr, freq, samplesInBaseFrequencyPeriod, TWO_PI);
	carrier(phasorFreqTable, sr, phasorFreq, samplesInPhasorFrequencyPeriod, TWO_PI);

	for(i = 0; i < end; i += blockframes)
	{
		for(j = 0; j < blockframes; j++, adsrIndex++, phaseIndex++)
		{
			if((phaseIndex + 1) >= samplesInBaseFrequencyPeriod)
			{
				phaseIndex = 0;
			}

			audioblock[j] = (invertedBaseFreqTable[phaseIndex % samplesInBaseFrequencyPeriod]) * 
				(phasorFreqTable[phaseIndex % samplesInPhasorFrequencyPeriod]) * 9000 * env[adsrIndex];

			// fills the rest of the file with 0's so the audio ends at
			// the end of the last wave that has completed the full wave/period/cycle
			if((adsrIndex + samplesInBaseFrequencyPeriod) > end)
			{
				for(; j < blockframes; j++)
				{
					audioblock[j] = 0;
				}
			}
		}
		// Append the created block to the .wav file
		// (Blocks are just an arbitrary number of samples, they
		// are not necessarily the start or the end of any wave)
		fwrite(audioblock, sizeof(short), blockframes, fpout);
	}

	// Prevent memory leaks and close the output .wav file
	// Memory leak prevention isn't strictly necessary, as
	// the program terminates right after and would thus free
	// the allocated memory, but it's a good habit none-the-less
	free(audioblock);
	free(header);
	fclose(fpout);

	return 0;
}

// This function creates the ADSR envelope. The array
// contains a value corresponding to the volume for the
// sample at a particular index. The index represents
// time that has elapsed, or in other words, the number
// of samples that have passed.
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

		for(int i = 0; i < dur * sr; i++)
		{
			arr[i] = 0.8;
		}
	}
	else
	{
		double increment = 1. / (period*a);
		for(int i = 0; i < period*a; i++)
		{
			arr[i] = increment*i;
		}

		increment = (1. - s) / (period*d);
		for(int i = period*a; i < period*(a+d); i++)
		{
			arr[i] = 1. - (increment*(i-(period*a)));
		}
		for(int i = period*(d+a); i < (sr*dur - period*r); i++)
		{
			arr[i] = s;
		}

		increment = s / (period*r);
		for(int i = period*r; i > 0; i--)
		{
			arr[sr*dur - i] = increment*i;
		}
	}
}

void modulator(double* invertedBaseFreqTable, int sr, int freq, int samplesInBaseFrequencyPeriod, double TWO_PI)
{
	for(int phaseIndex = 75, arrayPos = 0; phaseIndex < (samplesInBaseFrequencyPeriod + 75); phaseIndex++, arrayPos++)
	{
		double sample = 0;

		for(int n = 1; n < 200; n++)
		{
			if(n * freq >= sr / 2)	
			{ 
				// nyquist frequency check
				break;
			}
			
			// Sawtooth wave 
			sample += pow(-1, (n + 1)) / n * sin(phaseIndex * TWO_PI * n * freq / sr);
		}

		invertedBaseFreqTable[arrayPos] = (-1 * ((sample+1.5))) + 3.2;
	}
}

void carrier(double *arr, int sr, int phasorFreq, int samplesInPhasorFrequencyPeriod, double TWO_PI)
{
	for(int i = 0; i < samplesInPhasorFrequencyPeriod; i++)
	{
		arr[i] = sin(i*(TWO_PI/samplesInPhasorFrequencyPeriod));
	}
}

// This functon call creates the header for a .wav file
// It contains all the meta-data used by programs to decode
// the contents of the wave file. The audio data begins
// immediately after the end of the header.
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