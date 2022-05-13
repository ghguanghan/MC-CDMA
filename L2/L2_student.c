#include "OFDM.h"

/* Generates pilot sequence (random QPSK symbols)
 * Parameters
 *  pltI    - pilot sequence I-component samples
 *  pltQ    - pilot sequence Q-component samples
 *  len     - pilot sequence length in symbols
 */

void generate_pilot_qpsk(double pltI[], double pltQ[], int len) {
    // Randomly generate QPSK symbols

	//***Student_code_start*** 
    for(int i = 0; i < len; i++){
        double ra = (double) rand()/RAND_MAX; 
        pltI[i] = (ra > 0.5) ? SCALE_QPSK : -SCALE_QPSK;
        pltQ[i] = (ra > 0.5) ? SCALE_QPSK : -SCALE_QPSK;    
    }
	//***Student_code_end*****
}


/* Inserts TDM pilot sequence into the frame (OFDM symbol)
 * Parameters
 *  frmI    - frame sequence (OFDM symbol) I-component
 *  frmQ    - frame sequence (OFDM symbol) Q-component
 *  pltI    - pilot sequence I-component samples
 *  pltQ    - pilot sequence Q-component samples
 *  len     - pilot sequence length in symbols
 */

void insert_pilot_tdm(double frmI[], double frmQ[], double pltI[], double pltQ[], int len) {
    
    //insert_pilot_tdm(txSymI, txSymQ, pltI, pltQ, SYMBOLS_PER_BLOCK);
	//***Student_code_start*** 
    for(int i = 0; i < len; i++){
        frmI[i] = pltI[i];
        frmQ[i] = pltQ[i];
    }
	//***Student_code_end*****
}


/* Get data and pilot subcarrier indexes.
 *
 * Parameters
 *  dataIndx    - array to be filled with indexes of data subcarries
 *  pltIndx     - array to be filled with indexes of pilot subcarries
 *  numCarriers - total number of carriers (data + pilots)
 *  numPilots   - number of pilot subcarriers
 *  lastPilot   - indicator if the last carruer should be pilot (1) or not (0)
 *  
 * >>
 * NOTE: Edge pilots are better for polynomial interpolation as there is no
 *       extrapolation error at the edges. However, the high-resolution 
 *       interpolation based on FFT fails in that case.
 * >>
 * 
 */


void get_data_pilot_indexes(int dataIndx[], int pltIndx[], int numCarriers, int numPilots, int lastPilot){
	//***Student_code_start***
    int boundIdx = numCarriers - 1;
    int disP = boundIdx/(numPilots-1); //Calculation distance between two adjacent Pilots
    int disCounter = 0; // Distance Counter
    int dataCounter = 0;

    if(lastPilot == 1){ // last one must be pilot
        for(int i=0; i<numPilots; i++){
            /* Pilot Index */
            if(boundIdx % disP == 0){ //absolute equally spaced, disP = 3, 7, 9, 21
                pltIndx[i] = disCounter;
                disCounter += disP;
            }else{ // the Pilots are not absolute equally spaced
                if(i+1 == numPilots){
                    pltIndx[i] = boundIdx;
                }else{
                    pltIndx[i] = disCounter;
                    disCounter += disP;
                }
            }
            /* Data Index */
            while((dataCounter + i) != pltIndx[i] && i > 0){
                dataIndx[dataCounter] = i + dataCounter;
                dataCounter++;
            }
        }
    }else{ // last one is not pilot for FFT
        // number of carrier could be: 2, 4, 8, 16, 32
        int disPlt = numCarriers/numPilots;
        for(int i=1; i <= numPilots; i++){
            pltIndx[i-1] = disCounter;
            disCounter += disPlt;
            while(dataCounter + i < pltIndx[i-1] + disPlt){
                dataIndx[dataCounter] = i + dataCounter;
                dataCounter++;
            }
        }
    }
	//***Student_code_end*****lse dataIndx[j++] = i;

}


/* Insert symbols onto subcarriers with given indexes.
 *
 * Parameters
 *  frmI    - frame sequence (OFDM symbol) I-component
 *  frmQ    - frame sequence (OFDM symbol) Q-component
 *  symIndx - indexes of subcarries to insert symbols
 *  symI    - symbol sequence I-component samples
 *  symQ    - symbol sequence Q-component samples
 *  symLen  - symbol sequence length
 * 
 * >>
 * NOTE: Used for both pilot and data symbols.
 * >>
 */

void insert_symbols_fdm(double frmI[], double frmQ[], double symI[], double symQ[], int symIndx[], int symLen) {
	//***Student_code_start*** 
	for(int i = 0; i < symLen; i++){
		frmI[symIndx[i]] = symI[i];
		frmQ[symIndx[i]] = symQ[i];
	}
	//***Student_code_end*****
}

////////////////////
///Receiver/////////
////////////////////
//////////////////////////
//  CHANNEL ESTIMATION  //
//////////////////////////

/* Estimate channel (frequency) transfer function (CTF), for
 * TDM pilot, i.e. when all carriers in a frame are pilots.
 * 
 * Parameters
 *  Hi      - channel transfer function's I-component
 *  Hq      - channel transfer function's Q-component
 *  pltI    - Rx symbol sequence I-component
 *  pltQ    - Rx symbol sequence Q-component
 *  pltQ    - pilot sequence Q-component samples
 *  pltI    - pilot sequence I-component samples
 *  pltQ    - pilot sequence Q-component samples
 *  len     - pilot sequence length in symbols
 */


void ch_estimation(double Hi[], double Hq[], double rxSigI[], double rxSigQ[], double pltI[], double pltQ[], int len) {
	//***Student_code_start***
    // ch_estimation(HestI, HestQ, rxEstI, rxEstQ, pltI, pltQ, SYMBOLS_PER_BLOCK);
    for(int i = 0; i < len; i ++){
        double square_sum = pow(pltI[i],2) + pow(pltQ[i],2);
        Hi[i] = (rxSigI[i]*pltI[i] + rxSigQ[i]*pltQ[i])/square_sum;
        Hq[i] = (rxSigQ[i]*pltI[i] - rxSigI[i]*pltQ[i])/square_sum;
    }
	//***Student_code_end*****
}


/* Filter noise from the estimated CTF.
 *
 * Parameters
 *  Hi  - CTF I-component
 *  Hq  - CTF Q-component
 *  len - common sequence length
 * NOTE: Memory allocation on every run is not efficient.
 *       If implemented as a structure, preallocated array
 *       could be kept internally to avoid repeated allocations.
 */

void filter_noise(double Hi[], double Hq[], int len){
    // filter noise from CTF
    //***Student_code_start***
    double H_i[len];
    double H_q[len];
    for(int i = 0; i < len; i ++){    
        H_i[i] = 0;
        H_q[i] = 0;
        for( int j = 0; j < len; j++){
            double phase = 2*PI*j*i/len;
            H_i[i] += Hi[j]*cos(phase) - Hq[j]*sin(phase);
            H_q[i] += Hq[j]*cos(phase) + Hi[j]*sin(phase);
        }
        H_i[i] /= sqrt(len);
        H_q[i] /= sqrt(len);
        // Filter the Noise
        if(i >=  MAX_TAPS){
            H_i[i] = 0;
            H_q[i] = 0;
        }
    }
    fft(H_i, H_q, len, Hi, Hq);

	//***Student_code_end*****
}


// >> FDM PILOT <<

/* Estimate channel (frequency) transfer function (CTF)
 * at frequency points determined by pilot subcarrier indxes.
 *
 * Parameters
 *  Hi      - CTF I-component
 *  Hq      - CTF Q-component
 *  rxSigI  - Rx symbol sequence I-component
 *  rxSigQ  - Rx symbol sequence Q-component
 *  pltI    - pilot sequence I-component samples
 *  pltQ    - pilot sequence Q-component samples
 *  pltIndx - indexes of pilot subcarriers
 *  pltLen  - pilot sequence length in symbols
 * 
 * >>
 * NOTE: This function can be also used for TDM pilot, but an array
 *       with all indexes (0, 1, ..., numCarriers-1) has tu be provided
 *       in place of pltIndx.
 * >>
 */


void ch_estimation_fdm(double Hi[], double Hq[], double rxSigI[], double rxSigQ[], double pltI[], double pltQ[], int pltIndx[], int pltLen) {
	//***Student_code_start***
    // Channel Transfer Function (CTF) estimate
     for(int i = 0; i < pltLen; i++){
        double square_sum = pow(pltI[i],2) + pow(pltQ[i],2);
        Hi[pltIndx[i]] = (rxSigI[pltIndx[i]]*pltI[i] + rxSigQ[pltIndx[i]]*pltQ[i])/square_sum;
        Hq[pltIndx[i]] = (rxSigQ[pltIndx[i]]*pltI[i] - rxSigI[pltIndx[i]]*pltQ[i])/square_sum;
    }
	//***Student_code_end*****
}

/* Linear interpolation of the CTF based on available samples
 * at frequency points determined by pilot subcarrier indxes.
 *
 * Parameters
 *  Hi          - CTF I-component
 *  Hq          - CTF Q-component
 *  pltIndx     - indexes of pilot subcarriers
 *  numCarriers - total number of subcarriers (data + pilots)
 */


void interp_linear(double Hi[], double Hq[], int pltIndx[], int numCarriers, int numPilots) {
	//***Student_code_start*** 
    for(int numP = 1; numP < numPilots; numP++){
        int K = pltIndx[numP] - pltIndx[numP-1];
        for(int n = pltIndx[numP-1] + 1; n < pltIndx[numP]; n++){
             int i = n - K * (int) floor((double)n/K);
	     Hi[n] = (double) (K-i)/K * Hi[ pltIndx[ (int) floor( (double) n/K) ] ] + (double) i/K * Hi[ pltIndx[ (int) floor((double) n/K ) + 1] ];
	     Hq[n] = (double) (K-i)/K * Hq[ pltIndx[ (int) floor( (double) n/K) ] ] +  (double) i/K * Hq[ pltIndx[ (int) floor( (double) n/K ) + 1] ];
        }
    }/*
    int i, j, n, K;
	for(j=0; j<numPilots-1; j++){
        K = pltIndx[j+1] - pltIndx[j];
		for(n=pltIndx[j]+1; n<pltIndx[j+1]; n++){
            i = n - K * (int) floor((double)n/K);
			Hi[n] = (double) (K-i)/K * Hi[ pltIndx[ (int) floor( (double) n/K) ] ] + (double) i/K * Hi[ pltIndx[ (int) floor((double) n/K ) + 1] ];
			Hq[n] = (double) (K-i)/K * Hq[ pltIndx[ (int) floor( (double) n/K) ] ] +  (double) i/K * Hq[ pltIndx[ (int) floor( (double) n/K ) + 1] ];
		}
	}*/
	//***Student_code_end*****
}


/* Second-order polynomial (quadratic) interpolation of the CTF based on
 * available samples at frequency points determined by pilot subcarrier indxes.
 *
 * Parameters
 *  Hi          - CTF I-component
 *  Hq          - CTF Q-component
 *  pltIndx     - indexes of pilot subcarriers
 *  numCarriers - total number of subcarriers (data + pilots)
 */


void interp_quadratic(double Hi[], double Hq[], int pltIndx[], int numCarriers, int numPilots) {
  	//***Student_code_start*** 
    for(int numP = 1; numP < numPilots; numP++){
        int K = pltIndx[numP] - pltIndx[numP-1];
        for(int n = pltIndx[numP-1] + 1; n < pltIndx[numP]; n++){
            int i = n - K * (int) floor((double)n/K);
        Hi[n] = (double) (i - K)*(i - 2*K)/(2*K*K) * Hi[ pltIndx[ (int) floor( (double) n/K) ] ] +
                (double) i*(2*K - i)/(K*K) * Hi[ pltIndx[ (int) floor((double) n/K ) + 1] ] +
                (double) i*(i - K)/(2*K*K) * Hi[ pltIndx[ (int) floor((double) n/K ) + 2] ];
        Hq[n] = (double) (i - K)*(i - 2*K)/(2*K*K) * Hq[ pltIndx[ (int) floor( (double) n/K) ] ] +
                (double) i*(2*K - i)/(K*K) * Hq[ pltIndx[ (int) floor((double) n/K ) + 1] ] +
                (double) i*(i - K)/(2*K*K) * Hq[ pltIndx[ (int) floor((double) n/K ) + 2] ];
       }
   }
	//***Student_code_end*****
}


/* High-resolution interpolation based on FFT.
 *
 * Parameters
 *  Hi           - CTF sequence I-component samples
 *  Hq           - CTF sequence Q-component samples
 *  pltIndx      - indexes of pilot subcarriers
 *  numCarriers  - total number of subcarriers (data + pilots)
 *  numPilots    - total number of pilots
 */
// Interpolates the CTF at missing points by high-resolution DFT method.


void interp_fft(double Hi[], double Hq[], int pltIndx[], int numCarriers, int numPilots) {
  	//***Student_code_start*** 
	
    double hi[numCarriers];
    double hq[numCarriers];
    double hi_plt[numPilots];
    double hq_plt[numPilots];
    
    /* Step 1:
     * -> Using the index of the pilots, assign the pilots to those locations
     */
    for(int i=0; i < numPilots; i++){
        hi_plt[i] = Hi[pltIndx[i]];
        hq_plt[i] = Hq[pltIndx[i]];
    }
    
    /* Step 2:
     * -> transform the pilots points in to time domain
     * -> Filtering the carriers, which are larger than numPilots, to be ZERO
     */
    ifft(hi_plt, hq_plt, numPilots, hi, hq);
    int np = (numPilots < MAX_TAPS) ? numPilots : MAX_TAPS;
    for(int i = 0; i < numCarriers; i++){
        if(i >=  np){
            hi[i] = 0;
            hq[i] = 0;
        }
    }
    /* Step 3:
     * -> transform back to Frequency domain
     */
    fft(hi, hq, numCarriers, Hi, Hq);
	//***Student_code_end*****
	
}




