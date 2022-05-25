#include "OFDM.h"


//////////////////////////////////////
//  FORWARD ERROR CORRECTION (FEC)  //
//////////////////////////////////////

/* Encode 4-bit info block into an 8-bit codeword, by
 * applying a SYSTEMATIC Hamming (8,4) code.
 *
 * Parameters
 *  infoBits - input 4-bit sequence of information bits
 *  encBits  - output 8-bit sequence of encoded bits
 *
 */
void hamming84_encoder(unsigned char infoBits[], unsigned char encBits[]) {
    /*
     * Insert your code here
     */
    // 3567 0124
    //  v     v
    // 0123 4567
    encBits[5] = infoBits[0] ^ infoBits[1] ^ infoBits[3];
    encBits[6] = infoBits[0] ^ infoBits[2] ^ infoBits[3];
    encBits[7] = infoBits[1] ^ infoBits[2] ^ infoBits[3];
    encBits[4] = encBits[5] ^ encBits[6] ^ infoBits[0] ^ encBits[7] ^ infoBits[1] ^ infoBits[2] ^ infoBits[3];
    encBits[0] = infoBits[0];
    encBits[1] = infoBits[1];
    encBits[2] = infoBits[2];
    encBits[3] = infoBits[3];
}


/* Encode input bit sequence by applying Hamming (8,4) code over successive 4-bit blocks.
 *
 *  dataBits     - input (info) bit sequence
 *  encBits      - output (encoded) bit sequence
 *  infoBitsLen  - input sequence length (num. info bits)
 *
 */
void fec_encoder(unsigned char infoBits[], unsigned char encBits[], int infoBitsLen) {
    /*
     * Insert your code here
     */
    int numSym = infoBitsLen / QAM16; // 128/4 = 32 Symbol
    unsigned char tempInfoBits[QAM16];
    unsigned char tempEncBits[2*QAM16];
    for(int i = 0; i < numSym ; i++){//32
        for(int j = 0; j < QAM16; j ++){
            tempInfoBits[j] = infoBits[i*QAM16 + j];
        }
        hamming84_encoder(tempInfoBits, tempEncBits);
        for(int q = 0; q < 8; q++){
            encBits[i * 8 + q] = tempEncBits[q];
        }
    }
}


////////////////////
//  INTERLEAVING  //
////////////////////

/* Interleave the input sequence by writing cwLen bits into rows,
 * and reading out the obtain matrix column by column.
 *
 *  bitSeq     - input/output (interleaved) sequence
 *  mat     - interleaving matrix
 *  bitSeqLen     - bit sequence length
 *  numCols - codeword length of the applied FEC code
 *
 * >>
 * NOTE: Implementation with interleaving done in-place, i.e.
 *       result is placed back to the input array.
 * >>
 */
void bit_interleaver(unsigned char bitSeq[], unsigned char mat[], int bitSeqLen, int numCols) {
    /*
     * Insert your code here
     * => bit_interleaver(txBits, intMat, bitsPerBlock, cwLen);
     */
    int numRows = bitSeqLen/numCols;
    memcpy(mat, bitSeq, bitSeqLen*sizeof(unsigned char));
    for(int i = 0; i < numRows; i++){
        for(int j = 0; j < numCols; j++){
            bitSeq[i*numCols + j] = mat[j*numRows + i];
        }
    }
}

//////////////////////////////////////
//  FORWARD ERROR CORRECTION (FEC)  //
//////////////////////////////////////


/* Decode the SYSTEMATIC Hamming (8,4) code.
 *
 * Parameters
 *  infoBits - input 4-bit sequence of information bits
 *  encBits  - output 8-bit sequence of encoded bits
 *
 * Output
 *  int - indicator of number of errors, i.e.
 *          0 - no errors
 *          1 - single error (corrected)
 *          2 - double error (uncorrectable)
 *
 */
int hamming84_decoder(unsigned char encBits[], unsigned char infoBits[]) {
    /*
    * Insert your code here
    */
    unsigned char sydr[4];
    
    sydr[0] = encBits[4] ^ encBits[5] ^ encBits[6] ^ encBits[7] ^ encBits[0] ^ encBits[1] ^ encBits[2] ^ encBits[3];
    sydr[1] = encBits[5] ^ encBits[0] ^ encBits[1] ^ encBits[3];
    sydr[2] = encBits[6] ^ encBits[0] ^ encBits[2] ^ encBits[3];
    sydr[3] = encBits[7] ^ encBits[1] ^ encBits[2] ^ encBits[3];
    
    if(sydr[1] == 0 && sydr[2] == 0 && sydr[3] == 0){
        infoBits[0] = encBits[0];
        infoBits[1] = encBits[1];
        infoBits[2] = encBits[2];
        infoBits[3] = encBits[3];
        return 0; // No error
    }
    // 1 - single error
    if(sydr[1] && sydr[2] && sydr[3] == 0){
        infoBits[0] = !encBits[0];
        infoBits[1] = encBits[1];
        infoBits[2] = encBits[2];
        infoBits[3] = encBits[3];
        return 1;
    }else if(sydr[1] && sydr[3] && sydr[2] == 0){
        infoBits[0] = encBits[0];
        infoBits[1] = !encBits[1];
        infoBits[2] = encBits[2];
        infoBits[3] = encBits[3];
        return 1;
    }else if(sydr[2] && sydr[3] && sydr[1] == 0){
        infoBits[0] = encBits[0];
        infoBits[1] = encBits[1];
        infoBits[2] = !encBits[2];
        infoBits[3] = encBits[3];
        return 1;
    }else if(sydr[1] && sydr[2] && sydr[3]){
        infoBits[0] = encBits[0];
        infoBits[1] = encBits[1];
        infoBits[2] = encBits[2];
        infoBits[3] = !encBits[3];
        return 1;
    }else{
        // 2 Errors or 4 Errors
        infoBits[0] = encBits[0];
        infoBits[1] = encBits[1];
        infoBits[2] = encBits[2];
        infoBits[3] = encBits[3];
        return 2;
    }
}


/* Decode bit sequence encoded with Hamming (8,4) code,
 * and report if uncorrectable errors are detected.
 *
 * Parameters
 *  encBits  - input sequence of FEC encoded bits
 *  infoBits - output sequence of (decoded) info bits
 *  infoBitsLen  - output sequence length (num. info bits)
 *
 * Output
 *  ack_flg - successful transmission acknowledgment
 *              0 - no uncorrectable errors in the input sequence infoBits
 *              1 - at least one codeword is uncorrectable
 */
int fec_decoder(unsigned char encBits[], unsigned char infoBits[], int infoBitsLen) {
    /*
    * Insert your code here
    */
    // fec_decoder(rxBits, rxInfoBits, infoBitsPerBlock);
    int numSym = infoBitsLen / QAM16; // 128/4 = 32 Symbol
    int ct = 0;
    unsigned char tempInfoBits[QAM16];
    unsigned char tempEncBits[2*QAM16];
    for(int i = 0; i < numSym ; i++){//32
        for(int q = 0; q < 8; q++){
             tempEncBits[q] = encBits[i * 8 + q];
        }
        hamming84_decoder(tempEncBits, tempInfoBits);
        if(hamming84_decoder(tempEncBits, tempInfoBits) == 2){
            ct++;
        }
        for(int j = 0; j < QAM16; j ++){
            infoBits[i*QAM16 + j] = tempInfoBits[j];
        }
    }
    if(ct != 0) return 1;
    else return 0;
}



////////////////////
//  INTERLEAVING  //
////////////////////

/* Deinterleave the input sequence by reverting the interleaving operation.
 *
 *  bitSeq     - input/output (interleaved) sequence
 *  mat     - deinterleaving matrix
 *  bitSeqLen     - bit sequence length
 *  numCols - number of matrix columns or codeword length of the applied FEC code
 *
 * >>
 * NOTE: Deinterleaving is done in-place (result is saved to the same array).
 * >>
 */
void bit_deinterleaver(unsigned char bitSeq[], unsigned char mat[], int bitSeqLen, int numCols) {
    /*
     * Insert your code here
     * bit_deinterleaver(rxBits, deintMat, bitsPerBlock, cwLen);
     */
    int numRows = bitSeqLen/numCols;
    for(int i = 0; i < numRows; i++){
        for(int j = 0; j < numCols; j++){
            mat[j*numRows + i] = bitSeq[i*numCols + j];
        }
    }
    memcpy(bitSeq, mat, bitSeqLen*sizeof(unsigned char));
}

