/*
 COMPILE with:
    gcc -o L1_mccdma L1_mccdma.c Transmitter.c Receiver.c Channel.c Auxiliary.c -lm -Wall

 Usage example:
    ./L1_mccdma 2 5000 8

 Zip folder on server for download (call from 'Documents' folder):
    zip -r L1_mccdma.zip L1_mccdma
*/

/* Zheng
Compile with command line: 
	make

Usage example:
	./L1_mccdma 2 5000 8
*/

#include "OFDM.h"

void usage(char* progName) {
    printf("\nUsage: %s <taps> <symbol blocks> <spreading>\n\n", progName);
    printf("taps: length of the channel delay profile\n");
    printf("symbol blocks: number of symbol blocks to simulate\n");
    printf("spreading: [ 1 | 2 | 4 | ... | %d ]\n", SYMBOLS_PER_BLOCK);
}

void fix_file_name(char* fileName, char* modulation, char* howManyTaps,
                   char* spreading, char* kernel) {

    int lenMod = strlen(modulation);

    char tmp[FILE_NAME_SIZE];
    strcpy(tmp, kernel);
    strcat(tmp, "_");

    for (int i = 0; i < lenMod; i++) {
        modulation[i] = toupper( modulation[i] );
    }
    strcat(tmp, modulation);

    strcat(tmp, "_Taps");
    strcat(tmp, howManyTaps);
    strcat(tmp, "_SF");
    strcat(tmp, spreading);
    strcat(tmp, ".txt");

    strcpy(fileName, tmp);
}

//generate PN sequence

void pn_generator(char* pnSequence, int len) {

	//***Student_code_start***
    // Step 1: m = Number of Flip Flop(FF)
    int m = log2(len);
    char* tempPN = (char*)malloc(SYMBOLS_PER_BLOCK * m * sizeof(char));// there are m*symbol_per_block bits in tempPN array

    // Step 2: Inital all FF, as long as not: 000 000, it is ok.
    tempPN[0] = 1;
    tempPN[1] = 0;
    tempPN[2] = 0;
    tempPN[3] = 0;
    tempPN[4] = 0;
    tempPN[5] = 0;
    pnSequence[0] = tempPN[5];
    // Step 3: Build Logic, Start from second row.
    for(int i = 1; i < len; i++){
        tempPN[i*m] = tempPN[(i-1)*m] ^ tempPN[i*m-1];// calculate for 1st FF
        for(int j = 1; j < m; j++){
            tempPN[i*m + j] = tempPN[(i-1)*m + j-1]; // calculete shifted bits
        }
        // Step 4: Assign PN Sequence to actural pnSequence
        pnSequence[i] = tempPN[i*m - 1];
    }
    
	//***Student_code_end*****
    
}

int generate_symbols_mccdma(unsigned char* txBits, int bitsPerSymbol,
        int howManySymbols, int spreadFactor, double* txSymI, double* txSymQ){
        
	//***Student_code_start***
    /*
     * Modulation and Copier
     *
     */
    for(int i = 0; i < howManySymbols; i++){//Calculate I and Q Component
        // Encoder and Copier:
        for(int s = 0; s < spreadFactor; s++){
            txSymI[i*spreadFactor + s] = (txBits[i*bitsPerSymbol + 1]) ? SCALE_QPSK : -SCALE_QPSK;
            txSymQ[i*spreadFactor + s] = (txBits[i*bitsPerSymbol]) ? SCALE_QPSK : -SCALE_QPSK;
        }
    }
	//***Student_code_end*****

    return 0;
}

void spread_symbols(double* sigI, double* sigQ, char* pnSeq, int len) {

	//***Student_code_start*** 
	/*
     * For Spreading and Despreading
     *
     */
    for(int i = 0; i < len; i++){
        sigI[i] *= pnSeq[i];
        sigQ[i] *= pnSeq[i];
    }
    
	//***Student_code_end*****
}

void decode_symbols_mccdma(double* rxSymI, double* rxSymQ, int howManySymbols,
        int spreadFactor, int bitsPerSymbol,
        double sqrtSNR, unsigned char* rxBits) {

	//***Student_code_start***
    /*
     * For Decoding and Decopier and
     *
     */
    for(int i = 0; i < howManySymbols; i++){
        // Decopier
        double tempI = 0;
        double tempQ = 0;
        for(int s= 0; s < spreadFactor; s++){
            // Average all received copied I and Q samples
            tempI += rxSymI[i*spreadFactor + s]/(sqrtSNR * SQRT_OF_2*spreadFactor);
            tempQ += rxSymQ[i*spreadFactor + s]/(sqrtSNR * SQRT_OF_2*spreadFactor);
        }
        // Decoding
        rxBits[i*bitsPerSymbol+1] = (tempI > 0) ? 1 : 0;
        rxBits[i*bitsPerSymbol] = (tempQ > 0) ? 1 : 0;
    }
    
    
	//***Student_code_end*****


}

int main(int argc, char** argv) {

    /*
     * Initialize the random generator. Must be done once per run of main
     */
    srand((unsigned) time(NULL));

    double EbN0[] = {0, 5, 10, 15, 20, 25, 30};

    if (argc < 4) {
        usage(argv[0]);
        return(EXIT_FAILURE);
    }

    char modulation[] = "QPSK";
    int bitsPerSymbol = QPSK;

    int taps = atoi(argv[1]);

    int howManyBlocks = atoi(argv[2]);

    int spreadFactor = atoi(argv[3]);

    // ToDo: Check spread factor for correctness
    int howManySymbols = SYMBOLS_PER_BLOCK / spreadFactor;
    int howManyBits = howManySymbols * bitsPerSymbol;

    // Reuse to generate info bits at each iteration
    unsigned char txBits[howManyBits];
    unsigned char rxBits[howManyBits];

    /*
     * Losses due to GI insertion and spreading
     * For spreading it is simply spreadFactor
     */
    double cpLoss = (double) SYMBOLS_PER_BLOCK / (SYMBOLS_PER_BLOCK + CP_LEN);

    // SNR per bit
    double* snr = (double*) malloc(sizeof(EbN0));
    double* sqrtSNR = (double*) malloc(sizeof(EbN0));
    for (int i = 0; i < sizeof(EbN0) / sizeof(double); i++) {
        snr[i] = pow( 10, EbN0[i] / 10);

        /*
         * To normalize energy per bit. Models "modulation loss".
         * Makes comparison of different modulations fair.
         */
        snr[i] *= bitsPerSymbol;
        snr[i] /= spreadFactor;
        snr[i] *= cpLoss;
        sqrtSNR[i] = sqrt(snr[i]);
    }

    /*
     * We will generate tx symbols and trow them away in each iteration
     * Rx symbols will be stored for each block/iteration.
     * Add space for guard interval as well.
     */
    double* txSymI = (double*) malloc( SYMBOLS_PER_BLOCK * sizeof(double));
    double* txSymQ = (double*) malloc( SYMBOLS_PER_BLOCK * sizeof(double));

    /*
     * These are tx symbols after ifft is applied to modulate subcarriers
     */
    double* txModI = (double*) calloc( SYMBOLS_PER_BLOCK + CP_LEN, sizeof(double));
    double* txModQ = (double*) calloc( SYMBOLS_PER_BLOCK + CP_LEN, sizeof(double));

    double* rxSymI = (double*) calloc( SYMBOLS_PER_BLOCK + CP_LEN, sizeof(double));
    double* rxSymQ = (double*) calloc( SYMBOLS_PER_BLOCK + CP_LEN, sizeof(double));

    double* rxEstI = (double*) malloc((SYMBOLS_PER_BLOCK) * sizeof(double));
    double* rxEstQ = (double*) malloc((SYMBOLS_PER_BLOCK) * sizeof(double));

    /*
     * We will also generate channel taps for each iteration
     */
    double* hi = (double*) malloc(taps * sizeof(double));
    double* hq = (double*) malloc(taps * sizeof(double));

    double* alpha = (double*) malloc(N * taps * sizeof(double));
    double* phi = (double*) malloc(N * taps * sizeof(double));

    double* ber = (double*) malloc(sizeof(EbN0));

    /*
     * Spreading and despreading sequence
     * The length is equal to the number of subcarriers
     * Must be converted from [0,1] to [-1,1]
     */
    char* pnSeq = (char*)
            malloc(SYMBOLS_PER_BLOCK * sizeof(char));

    pn_generator(pnSeq, SYMBOLS_PER_BLOCK);

    for (int i = 0; i < SYMBOLS_PER_BLOCK; i++) {
        pnSeq[i] = 2 * pnSeq[i] - 1;
    }

    // >> SIMULATION <<

    // terminal output header
    printf("\nEbNo\tBER\n");

    // Repeat simulation for each given SNR point
    for (int i = 0; i < sizeof(EbN0) / sizeof(double); i++) {

        // initialize bit/block error counters
        ber[i] = 0;


        for (int j = 0; j < howManyBlocks; j++) {

            ///////////////////
            //  TRANSMITTER  //
            ///////////////////


            // generate data bits
            generate_info_bits(howManyBits, txBits);

            // modulate bit sequence
            generate_symbols_mccdma(txBits, bitsPerSymbol, howManySymbols,
                    spreadFactor, txSymI, txSymQ);

            // spreading symbols
            spread_symbols(txSymI, txSymQ, pnSeq, SYMBOLS_PER_BLOCK);

            // OFDM modulation
            ifft(txSymI, txSymQ, SYMBOLS_PER_BLOCK,
                    txModI + CP_LEN, txModQ + CP_LEN);

            // insert cyclic prefix
            insert_cp(SYMBOLS_PER_BLOCK, CP_LEN, txModI, txModQ);

            // scale Tx symbols to match SNR
            path_loss(txModI, txModQ, SYMBOLS_PER_BLOCK + CP_LEN, sqrtSNR[i]);

            ///////////////
            //  CHANNEL  //
            ///////////////

            // Fading
            multipath_fading_init(N, taps, fDT, alpha, phi);
            // We do not need correlated fading so 'time' must be fixed
            multipath_fading(alpha, phi, N, taps, 0, hi, hq);

            // Convolution with multipath
            conv_with_gi(hi, hq, taps, txModI, txModQ, SYMBOLS_PER_BLOCK + CP_LEN,
                    CP_LEN, rxSymI, rxSymQ);

            // AWGN
            add_noise(rxSymI, rxSymQ, SYMBOLS_PER_BLOCK + CP_LEN);

            ////////////////
            //  RECEIVER  //
            ////////////////

            /*
             * We will drop the GI by simply shifting by CP_LEN
             * and demodulate the subcarriers using fft
             */
            // OFDM demodulation (CI is dropped)
            fft(rxSymI + CP_LEN, rxSymQ + CP_LEN, SYMBOLS_PER_BLOCK, rxEstI, rxEstQ);

            // Equalization
            equalizerType equalizer = MMSE;
            //fde(equalizer, snr[i], hi, hq, rxEstI, rxEstQ, SYMBOLS_PER_BLOCK);
            fde(equalizer, snr[i], hi, hq, taps, rxEstI, rxEstQ, SYMBOLS_PER_BLOCK);

            // despreading symbols
            spread_symbols(rxEstI, rxEstQ, pnSeq, SYMBOLS_PER_BLOCK);

            // Demodulate (detect symbols and output data bits)
            decode_symbols_mccdma(rxEstI, rxEstQ, howManySymbols, spreadFactor,
                    bitsPerSymbol, sqrtSNR[i], rxBits);

            // Count bit and symbol errors in the block (only for data symbols)
            for (int m = 0; m < howManyBits; m++) {
                if (txBits[m] != rxBits[m])
                    ber[i]++;
            }
        }

        // Calculate BER
        ber[i] /= (howManyBits * howManyBlocks);
        printf("%f\t%f\n", EbN0[i], ber[i]);
    }

    // generate output file name
    char fileNameBER[FILE_NAME_SIZE];
    fix_file_name(fileNameBER, modulation, argv[1], argv[3], FILE_BER);

    // save BER to file
    save_asignal(EbN0, ber, sizeof(EbN0)/sizeof(double), fileNameBER);
    printf("\nBER saved to %s\n", fileNameBER);

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
    free(pnSeq);
    free(PN);
    
    return EXIT_SUCCESS;
}

