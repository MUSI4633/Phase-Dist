///Project 3: Audio Software Development
///Keegan Herperger 001180181
///Program to generate a tone in one of three wavehapes, with ADSR amplitude envelope.
///
///I couldn't get bit depth working; it changes the length of the file.
///A bit-depth of 16 will give the correct lengtht though.

#include <math.h>
#include <stdio.h> 
#include <stdlib.h>
#include <string.h>
#include "wave.h"  /* header file containing the wavehead struc 
					  and the declaration of update_header() */
FILE* fpout;

void modgen(int end, int blockframes, short* audioblock, int freq, int sr, int dur, double* env, int bd);
void breaks(float* plist, float* llist);
void addsyn(int dur);
void envelope(double *arr, int sr, int dur);
void carriergen(int end, int blockframes, short* audioblock, int freq, int sr, int dur, double* env, int bd);

int main(int argc, char** argv)
{
	short  *audioblock;      /* audio memory pointer */
	int     end, i, j;       /* dur in frames, counter vars */
	int     sr = 44100;      /* sampling rate */
	int     blockframes = 256; /* audio block size in frames */
	int     databytes;  /* audio data in bytes */
	unsigned int ndx = 0;   /* phase index for synthesis */
	float   dur, freq, sustain; /* duration, frequency, sustain */
	int attack, decay, release, depth, shape;   /* attack, decay, release, bit-depth, shape*/
	double  twopi;          /* 2*PI  */
	wavehead *header;

	if(argc != 6)
	{
		printf("usage: %s outfile dur freq bit-depth samplerate\n", argv[0]);
		exit(-1);
	}

	dur = atof(argv[2]);
	freq = atof(argv[3]);
	twopi = 8 * atan(1.);
	depth = atoi(argv[4]);
	sr = atof(argv[5]);
	end = (int)(dur*sr);
	fpout = fopen(argv[1], "wb");
	audioblock = (short *)malloc(sizeof(short)*blockframes);

	/* set the data size */
	databytes = end * sizeof(short);

	//allocate space for the envelope
	double *env;
	env = malloc(end * sizeof(double));

	/* write the header */
	header = (wavehead *)malloc(sizeof(wavehead));
	update_header(header, sr, 1, 16, databytes);
	fwrite(header, 1, sizeof(wavehead), fpout);

	// #TODO call implementation here
	
	free(audioblock);
	free(header);
	fclose(fpout);

	return 0;
}

void modgen(int end, int blockframes, short* audioblock, int freq, int sr, int dur, double* env, int bd)
{
	int jnum;
	printf("How many jumps would you like in the modulator signal?");
	scanf("%d", jnum);

	float* llist[jnum], plist[jnum];

	breaks(plist, llist); //we need to calculate, based on samplerate, how steep the
						  //linear segments need to be
	int sps = sr*(1 / freq);		//samples per cycle

	double twopi = 8 * atan(1.);
	int i, j, n, ndx = 0;

	envelope(env, sr, dur);

	// experimenting with new write function here
	for(i = 0; i < end; i += blockframes)
	{
		for(j = 0; j < blockframes; j++, ndx++)
		{
			audioblock[j] = 0;//8000*sin(ndx*twopi*freq/sr);

			for(int h = 1; h < 200; h++)
			{ //this is adding sinudoidal harmonics
				if(h*freq >= sr / 2)
				{ //nyquist frequency check
					break;
				}

				audioblock[j] += 16000 * (1. / h)*sin((ndx*twopi*h*freq) / sr);
			}
			if(ndx >= dur*sr)
			{ //this is ending the loop before it goes out of bounds
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

void breaks(float* plist, float* llist)
//int end, int blockframes, 
//		short* audioblock, int sr, int dur)
{
	bool flag = 0;
	int freq;

	printf("Input phasor frequency:");
	scanf("%d", freq);

	for(int i < sizeof(&plist); i = 0; i++)
	{
		printf("Input next jump point as a percentage of phase (0.-1.):");
		scanf("%s", &plist);

		printf("Input jump destination as a percentage of amplitude (0.-1.):");
		scanf("%s", &llist);

	}

}

void addsyn(int dur)
{
	int harmnum, relamp;

	printf("How many harmonics would you like in the carrier?");
	scanf("%d", &harmnum);

	printf("Choose the relative harmonic strength:\n 1. 1/n\n 2. 1/n^2\n 3.sqrt(n) (normalised)\n");
	scanf("%d", &relamp);

	carriergen(dur, harmnum, relamp);
}

// Function to generate the amplitude envelope.
//will give an error if it is too long.
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

//function to generate the square wave.
//will generate up to 200 harmonics, or until
//it reaches the Nyquist frequency.
void carriergen(int end, int blockframes, short* audioblock, int freq, int sr, int dur, double* env, int bd)
{
	double twopi = 8 * atan(1.);
	int i, j, n, ndx = 0;
	envelope(env, sr, dur);

	for(i = 0; i < end; i += blockframes)
	{
		for(j = 0; j < blockframes; j++, ndx++)
		{
			audioblock[j] = 0;//8000*sin(ndx*twopi*freq/sr);

			for(int h = 1; h < 200; h++)
			{
				n = (2 * h) - 1;
				if(n*freq >= sr / 2)
				{
					break;
				}

				audioblock[j] += 16000 * (1. / n)*sin((ndx*twopi*n*freq) / sr);
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