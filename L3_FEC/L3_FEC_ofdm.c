/*
 COMPILE with:
    gcc -o L3_FEC_ofdm L3_FEC_ofdm.c Transmitter.c Receiver.c Channel.c Auxiliary.c -lm -Wall

 Usage example:
    ./L3_FEC_ofdm 1 20000 1

*/

/* Kopic
Compile with command line: 
	make

Usage example:
	./L3_FEC_ofdm 1 20000 1
*/

#include "OFDM.h"


void usage(char* progName) {
    printf("\nUsage: %s <taps> <symbol blocks> <fec flag>\n\n", progName);
    printf("taps: length of the channel delay profile\n");
    printf("symbol blocks: number of symbol blocks to simulate\n");
    printf("FEC activation flag:\n");
    printf("   1 - Hamming (8,4) applied;\n");
    printf("   0 - No FEC applied.\n\n");
}

void output_file_name(char const* filePath, char* fileName, 
        char const* numTaps, int FECflg) {

    strcpy(fileName, filePath);
    strcat(fileName, "BER_QAM16_MMSE");

    strcat(fileName, "_L");
    strcat(fileName, numTaps);


    if (FECflg) {
        strcat(fileName, "_fec.txt");
    }
    else {
        strcat(fileName, "_nofec.txt");
    }
}

int main(int argc, char** argv) {
    
    // Initialize the random generator. Must be done once per run of main
    srand((unsigned) time(NULL));
    
    // set if SNR-per-bit values for BER/SER evalation
    double EbN0[] = {0, 5, 10, 15, 20, 25, 30, 35, 40};

    // check required if all required parameters are provided
    if (argc != 4) {
        usage(argv[0]);
        return(EXIT_FAILURE);
    }
    
    // set number of bits per symbol (corresponding to QPSK)
    int bitsPerSymbol = QAM16;

    // set equalizer type (MMSE)
    equalizerType equalizer = MMSE;
    
    // read the number of channel taps from parameter list
    int numTaps = atoi(argv[1]);

    // set number of blocks to simulate
    int numBlocks = atoi(argv[2]);

    // set FEC flag
    int FECflg = atoi(argv[3]);


    // Hammning (8,4) code parameters >
    int cwLen = 8; // codeword length
    double codeRate = 0.5; //  4/8

    // number of bits per block (OFDM symbol)
    int bitsPerBlock = SYMBOLS_PER_BLOCK * bitsPerSymbol; // total (info+redundancy) 
    
    // number of information bits
    int infoBitsPerBlock = bitsPerBlock;
    if (FECflg) infoBitsPerBlock *= codeRate; // considers code redundancy


    // loss due to GI insertion
    double cpLoss = (double) SYMBOLS_PER_BLOCK / (SYMBOLS_PER_BLOCK + CP_LEN);

    // loss due to FEC redundancy
    double fecLoss = 1.0;           // uncoded (energy used only for info bits)
    if (FECflg) fecLoss = codeRate; // coded (energy goes to redundancy bits as well)
 
    // Signal magnitude adjustment to match specified EbNo
    // NOTE: Noise is assumed to be of unit power per I/Q component
    double* snr = (double*) malloc(sizeof(EbN0));
    double* sqrtSNR = (double*) malloc(sizeof(EbN0));
    for (int i = 0; i < sizeof(EbN0) / sizeof(double); i++) {
        snr[i] = pow( 10, EbN0[i] / 10); // SNR per bit
        
        // Adjust SNR per symbol to match Eb/No
        snr[i] *= bitsPerSymbol;    // SNR per symbol
        snr[i] *= cpLoss;           // GI (cyclic prefix) insertion
        snr[i] *= fecLoss;          // FEC (redundancy insertion)
        sqrtSNR[i] = sqrt(snr[i]);  // amplitude scaling
    }

    // arrays to keep Tx and Rx information bits in each iteration (i.e. one OFDM symbol)
    unsigned char* txInfoBits = (unsigned char*) malloc( infoBitsPerBlock * sizeof(unsigned char));
    unsigned char* rxInfoBits = (unsigned char*) malloc( infoBitsPerBlock * sizeof(unsigned char));

    // encoded bits
    unsigned char* txBits = (unsigned char*) malloc( bitsPerBlock * sizeof(unsigned char));
    unsigned char* rxBits = (unsigned char*) malloc( bitsPerBlock * sizeof(unsigned char));

    // interleaver and deinterleaver matrix
    unsigned char* intMat = (unsigned char*) malloc( bitsPerBlock * sizeof(unsigned char));
    unsigned char* deintMat = (unsigned char*) malloc( bitsPerBlock * sizeof(unsigned char));

    // Tx symbols
    double* txSymI = (double*) malloc( SYMBOLS_PER_BLOCK * sizeof(double));
    double* txSymQ = (double*) malloc( SYMBOLS_PER_BLOCK * sizeof(double));
    
    // Tx OFDM symbol (modulated)
    double* txModI = (double*) calloc( SYMBOLS_PER_BLOCK + CP_LEN, sizeof(double));
    double* txModQ = (double*) calloc( SYMBOLS_PER_BLOCK + CP_LEN, sizeof(double));
    
    // Rx OFDM symbol (after channel)
    double* rxSymI = (double*) calloc( SYMBOLS_PER_BLOCK + CP_LEN, sizeof(double));
    double* rxSymQ = (double*) calloc( SYMBOLS_PER_BLOCK + CP_LEN, sizeof(double));
    
    // Rx symbols
    double* rxEstI = (double*) malloc((SYMBOLS_PER_BLOCK) * sizeof(double));
    double* rxEstQ = (double*) malloc((SYMBOLS_PER_BLOCK) * sizeof(double));
    
    // multipath channel
    double* hi = (double*) malloc(numTaps * sizeof(double));
    double* hq = (double*) malloc(numTaps * sizeof(double));
    
    double* alpha = (double*) malloc(N * numTaps * sizeof(double));
    double* phi = (double*) malloc(N * numTaps * sizeof(double));

    // Channel transfer function estimate
    double* HestI = (double*) malloc((SYMBOLS_PER_BLOCK) * sizeof(double));
    double* HestQ = (double*) malloc((SYMBOLS_PER_BLOCK) * sizeof(double));

    
    // Bit and symbol error calculation
    double* ber = (double*) malloc(sizeof(EbN0));


    // >> SIMULATION <<

    // terminal output header
    printf("\nEbNo\tBER\n");

    // Repeat simulation for each given SNR point
    for (int i = 0; i < sizeof(EbN0) / sizeof(double); i++) {
        
        // initialize bit/block error counters
        ber[i] = 0.0;

        // Initialize/reset multipath fading simulator
        multipath_fading_init(N, numTaps, fDT, alpha, phi);

        // OFDM simulation loop (block by block)
        for (int j = 0; j < numBlocks; j++) {

            ///////////////////
            //  TRANSMITTER  //
            ///////////////////

            // generate information bits
            generate_info_bits(infoBitsPerBlock, txInfoBits);

            if (FECflg) {
                // FEC encoder
                fec_encoder(txInfoBits, txBits, infoBitsPerBlock);
                
                // interleaving
                bit_interleaver(txBits, intMat, bitsPerBlock, cwLen);
            }
            else {
                // uncoded case
                memcpy(txBits, txInfoBits, infoBitsPerBlock*sizeof(unsigned char));
            }
            
            // modulate bit sequence
            generate_symbols(txBits, bitsPerSymbol, SYMBOLS_PER_BLOCK, txSymI, txSymQ);

            // OFDM modType
            ifft(txSymI, txSymQ, SYMBOLS_PER_BLOCK, txModI + CP_LEN, txModQ + CP_LEN);
            
            // insert cyclic prefix
            insert_cp(SYMBOLS_PER_BLOCK, CP_LEN, txModI, txModQ);
            

            // scale Tx symbols to match SNR                    
            path_loss(txModI, txModQ, SYMBOLS_PER_BLOCK + CP_LEN, sqrtSNR[i]);
            
            ///////////////
            //  CHANNEL  //
            ///////////////

            // MULTIPATH FADING
            // update channel impulse response at the beginning of each block
            multipath_fading(alpha, phi, N, numTaps, j*(SYMBOLS_PER_BLOCK + CP_LEN), hi, hq);
            
            // Convolution with channel impulse response
            conv_with_gi(hi, hq, numTaps, txModI, txModQ, SYMBOLS_PER_BLOCK + CP_LEN, 
                    CP_LEN, rxSymI, rxSymQ);
            
            // AWGN
            add_noise(rxSymI, rxSymQ, SYMBOLS_PER_BLOCK + CP_LEN);
            
            ////////////////
            //  RECEIVER  //
            ////////////////

            // Automatic Gain Control (AGC)
            // SNR estimation is typically done based on preamble (perfect estimate is assumed here)
            // NOTE: Reverts the pilot scaling introduced in path_loss function
            for (int k = CP_LEN; k < SYMBOLS_PER_BLOCK+CP_LEN; k++) {
                rxSymI[k] /= (sqrtSNR[i]*SQRT_OF_2);
                rxSymQ[k] /= (sqrtSNR[i]*SQRT_OF_2);
            }

            // OFDM demodulation (CI is dropped)
            fft(rxSymI + CP_LEN, rxSymQ + CP_LEN, SYMBOLS_PER_BLOCK, rxEstI, rxEstQ);

            // CTF for assumed perfect channel knowledge
            for (int k = 0; k < SYMBOLS_PER_BLOCK; k++) {
                HestI[k] = 0;
                HestQ[k] = 0;
                for (int m = 0; m < numTaps; m++) {
                    HestI[k] += hi[m] * cos(-2 * PI * k * m / SYMBOLS_PER_BLOCK) 
                            - hq[m] * sin(-2 * PI * k * m / SYMBOLS_PER_BLOCK);
                    HestQ[k] += hq[m] * cos(-2 * PI * k * m / SYMBOLS_PER_BLOCK) 
                            + hi[m] * sin(-2 * PI * k * m / SYMBOLS_PER_BLOCK);
                }
            }

            // Equalization (takes channel transfer function as input)
            fde(equalizer, snr[i], HestI, HestQ, rxEstI, rxEstQ, SYMBOLS_PER_BLOCK);
            
            // Demodulate (detect symbols and output data bits)
            decode_symbols(rxEstI, rxEstQ, SYMBOLS_PER_BLOCK, bitsPerSymbol, rxBits);

            if (FECflg) {
                // deinterleaving
                bit_deinterleaver(rxBits, deintMat, bitsPerBlock, cwLen);

                // FEC decoder
                fec_decoder(rxBits, rxInfoBits, infoBitsPerBlock);
            }
            else {
                // uncoded case
                memcpy(rxInfoBits, rxBits, infoBitsPerBlock*sizeof(unsigned char));
            }
            
            // count bit errors in the block
            for (int m = 0; m < infoBitsPerBlock; m++) {
                if (rxInfoBits[m] != txInfoBits[m]) ber[i]++;
            }
        }
            
        // Calculate BER
        ber[i] /= (infoBitsPerBlock * numBlocks);
        printf("%.1f\t%f\n", EbN0[i], ber[i]);
    }
    
    // generate output file name
    char fileNameBER[FILE_NAME_SIZE];
    output_file_name("./results/", fileNameBER, argv[1], FECflg);
    
    // save BER to file
    save_asignal(EbN0, ber, sizeof(EbN0)/sizeof(double), fileNameBER);
    printf("\nBER saved to %s\n", fileNameBER);
    
    // Release all allocated memory
    free(txInfoBits);
    free(rxInfoBits);
    free(txBits);
    free(rxBits);
    free(intMat);
    free(deintMat);
    free(txSymI);
    free(txSymQ);
    free(txModI);
    free(txModQ);
    free(rxSymI);
    free(rxSymQ);
    free(rxEstI);
    free(rxEstQ);
    free(snr);
    free(sqrtSNR);
    free(alpha);
    free(phi);
    free(hi);
    free(hq);
    free(ber);
    //
    free(HestI);
    free(HestQ);
    //
    return EXIT_SUCCESS;
}

