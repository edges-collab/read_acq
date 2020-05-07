/*=================================================================
 *
 * decode.c     Decodes the EDGES raw spectrum into double.
 *
 *
 * This is a MEX-file for MATLAB.
 *
 *=================================================================*/


#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

unsigned int decode_table[256];
char encode_table[64];

unsigned short decode_lookup_inited=0;
unsigned short encode_lookup_inited=0;

/* This makes a lookup table to convert an encoded char back into a
6-bit value */
static void initDecodeTable()
{
    int i = 0;
    for (i = 0; i<256; i++)
    {
        decode_table[i] = 0;
        if (i >= 'A' && i <= 'Z') {
            decode_table[i] = i - 'A';
        }
        else if (i >= 'a' && i <= 'z') {
            decode_table[i] = i - 'a' + 26;
        }
        else if (i >= '0' && i <= '9') {
            decode_table[i] = i - '0' + 52;
        }
        else if (i == '+') {
            decode_table[i] = 62;
        }
        else if (i == '/') {
            decode_table[i] = 63;
        }
    }
    decode_lookup_inited=1;
}

static void initEncodeTable()
{
    int i;

    // Define the encoding lookup table
    for (i = 0; i < 26; i++) {
        encode_table[i] = 'A' + i;
        encode_table[i + 26] = 'a' + i;
    }

    for (i = 0; i < 10; i++) {
        encode_table[i + 52] = '0' + i;
    }

    encode_table[62] = '+';
    encode_table[63] = '/';

    encode_lookup_inited=1;
}


int get_decoded_length(char *contents){
    return(strlen(contents)/4);
}

int decode(char *contents, double *out){
    /*
        Parameters
        ----------
        contents : a string of binary data from an ACQ file. Should be a single line
                   of the datafile, not the entire file.
        out : the decoded data contents. It should be pre-allocated before being
              sent to this function. Use the get_decoded_length function to
              determine the relevant length of out.
    */

    int decoded_length = get_decoded_length(contents);
    int i=0;

    /* Initialize the lookup table */
    if(!decode_lookup_inited){
        initDecodeTable();
    }

    /* Decode the data */
    for (i=0; i<decoded_length; i++)
    {
        int m = 4 * i;
        unsigned int temp = (decode_table[contents[m]] << 18) +
                            (decode_table[contents[m+1]] << 12) +
                            (decode_table[contents[m+2]] << 6) +
                            (decode_table[contents[m+3]]);

        /* Put the spectrum into linear units (it was stored in dB*-1e5) */
        out[i] = pow((double) 10.0, (double) -0.1 * (double) temp * (double) 0.00001);

    }

    return(0);
}

int encode(int length, double *in, char* out)
{
    int i, j, k;
    double dOut;
    int cnt= 0;

    if(!encode_lookup_inited) initEncodeTable();

    for(i=0; i<length; i++){

        // Convert power to decibels, or if power is zero (or at the start of the
        // spectrum), set to a negative number.
        if (in[i] <= 0) {
            dOut = -199;
        } else {
            dOut = 10 * log10(in[i]);
        }

        // Scale the decibel power to improve precision.
        k = -(int)(dOut * 100000);

        // Clip the power to relevant base64 range.
        if (k > 16700000) k = 16700000;
        else if (k < 0) k = 0;

        // Perform the encoding.
        for (j=0; j<4; j++) {
            out[cnt] = encode_table[(k >> (18-j*6)) & 0x3f];
            cnt += 1;
        }
    }

    return(0);

}
