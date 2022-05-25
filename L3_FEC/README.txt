1. Implementations of new functions go to "L3_Student.c", as indicated during the class.

2. Interfaces of new functions to be implemented for this programing exercise are placed in "OFDM.h".

    // FORWARD ERROR CORRECTION
    // Tx side
    void hamming84_encoder(unsigned char infoBits[], unsigned char encBits[]);
    void fec_encoder(unsigned char infoBits[], unsigned char encBits[], int infoBitsLen);

    // Rx side
    int hamming84_decoder(unsigned char encBits[], unsigned char infoBits[]);
    int fec_decoder(unsigned char encBits[], unsigned char infoBits[], int infoBitsLen);

    // INTERLEAVING
    // Tx side
    void bit_interleaver(unsigned char bitSeq[], unsigned char mat[], int bitSeqLen, int cwLen);

    // Rx side
    void bit_deinterleaver(unsigned char bitSeq[], unsigned char mat[], int bitSeqLen, int cwLen);


3. Please make sure you have the latest version of "gcc"

4. Compile with command line:
	make

5. Execute with command line.
	./L3_FEC_ofdm 8 20000 1
   E.g.: 8 means taps - the length of the channel delay profile
         20000 means symbol blocks - the number of symbol blocks to simulate
         1 means FEC activation flag:
            1 - Hamming (8,4) applied
            0 - No FEC applied
   Your results will be stored in a .txt file in 'results' folder.

6. Plot your graphs with any tools you like (e.g., Python, ...)
