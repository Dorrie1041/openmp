#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <iostream>

void SmithWaterman(int8_t *seq1,
                   int8_t *seq2,
                   const int32_t &length1,
                   const int32_t &length2,
                   const int &_match_score,
                   const int &_mismatch_score,
                   const int &_open_gap_penalty,
                   const int &_extend_gap_penalty,
                   int32_t &maxScore,
                   int32_t &endPosition1,
                   int32_t &endPosition2)
{

    int32_t **M = new int32_t *[length1 + 1];
    int32_t **E = new int32_t *[length1 + 1];
    int32_t **F = new int32_t *[length1 + 1];

    int32_t i, j;
    for (i = 0; i < (length1 + 1); ++i)
    {
        M[i] = new int32_t[length2 + 1];
        E[i] = new int32_t[length2 + 1];
        F[i] = new int32_t[length2 + 1];
        std::fill_n(M[i], length2 + 1, 0);
        std::fill_n(E[i], length2 + 1, 0);
        std::fill_n(F[i], length2 + 1, 0);
    }

    endPosition1 = 0;
    endPosition2 = 0;
    maxScore = 0;
    for (i = 1; i <= length1; ++i)
    {
        for (j = 1; j <= length2; ++j)
        {
            E[i][j] = (_open_gap_penalty + M[i][j - 1]) > (_extend_gap_penalty + E[i][j - 1]) ? (_open_gap_penalty + M[i][j - 1]) : (_extend_gap_penalty + E[i][j - 1]);
            //            E[i][j] = E[i][j] > (_open_gap_penalty + F[i][j-1]) ? E[i][j] : (_open_gap_penalty + F[i][j-1]);

            F[i][j] = (_open_gap_penalty + M[i - 1][j]) > (_extend_gap_penalty + F[i - 1][j]) ? (_open_gap_penalty + M[i - 1][j]) : (_extend_gap_penalty + F[i - 1][j]);
            //            F[i][j] = F[i][j] > (_open_gap_penalty + E[i-1][j]) ? F[i][j] : (_open_gap_penalty + E[i-1][j]);

            int32_t match = (seq1[i - 1] == seq2[j - 1]) ? (_match_score) : (_mismatch_score);
            match = match + M[i - 1][j - 1];

            M[i][j] = (match > E[i][j]) ? (match) : E[i][j];
            M[i][j] = M[i][j] > F[i][j] ? M[i][j] : F[i][j];
            M[i][j] = M[i][j] > 0 ? M[i][j] : 0;
            F[i][j] = F[i][j] > 0 ? F[i][j] : 0;
            E[i][j] = E[i][j] > 0 ? E[i][j] : 0;

            if (M[i][j] > maxScore)
            { // please do not change > to >=, since we are doing local alignment
                // >= will omit the first similar fragments
                maxScore = M[i][j];
                endPosition1 = i; // this means the endPosition1 and endPosition2 is 1 based coordinate
                endPosition2 = j;
            }
        }
    }
}

void SmithWaterman_Diag(int8_t *seq1,
                        int8_t *seq2,
                        const int32_t &length1,
                        const int32_t &length2,
                        const int &_match_score,
                        const int &_mismatch_score,
                        const int &_open_gap_penalty,
                        const int &_extend_gap_penalty,
                        int32_t &maxScore,
                        int32_t &endPosition1,
                        int32_t &endPosition2){

    int32_t **M = new int32_t *[length1 + 1];
    int32_t **E = new int32_t *[length1 + 1];
    int32_t **F = new int32_t *[length1 + 1];

    int32_t i, j;
    for (i = 0; i < (length1 + 1); ++i)
    {
        M[i] = new int32_t[length2 + 1];
        E[i] = new int32_t[length2 + 1];
        F[i] = new int32_t[length2 + 1];
        std::fill_n(M[i], length2 + 1, 0);
        std::fill_n(E[i], length2 + 1, 0);
        std::fill_n(F[i], length2 + 1, 0);
    }

    endPosition1 = 0;
    endPosition2 = 0;
    maxScore = 0;

    int32_t intervals = 1000;
    int32_t rows_chunks = (length1 + intervals -1) /intervals;
    int32_t cols_chunks = (length2 + intervals -1) /intervals;
    int32_t a, b;
    for (int32_t diag = 0; diag < rows_chunks + cols_chunks -1; ++diag){
         if (diag < cols_chunks){
            a = 0;
            b = diag;
         } else {
            a = diag - cols_chunks + 1;
            b = cols_chunks - 1; 
         }
         //do parallel here
         while (a < rows_chunks && b >= 0){
            for (i = a * intervals; i <= std::min((a + 1) * intervals, length1); ++i){
                for (j = b * intervals; j <= std::min((b + 1) * intervals, length2); ++j){
                        if (i == 0 || j == 0) continue; // Skip processing for i=0 or j=0
                        E[i][j] = (_open_gap_penalty + M[i][j-1]) > (_extend_gap_penalty + E[i][j-1]) ? (_open_gap_penalty + M[i][j-1]) : (_extend_gap_penalty + E[i][j-1]);
                        F[i][j] = (_open_gap_penalty + M[i-1][j]) > (_extend_gap_penalty + F[i-1][j]) ?( _open_gap_penalty + M[i-1][j]) : (_extend_gap_penalty + F[i-1][j]);
                        
                        int32_t match = (seq1[i - 1] == seq2[j - 1]) ? (_match_score) : (_mismatch_score);
                        match = match + M[i - 1][j - 1];

                        M[i][j] = (match > E[i][j]) ? (match) : E[i][j];
                        M[i][j] = M[i][j] > F[i][j] ? M[i][j] : F[i][j];
                        M[i][j] = M[i][j] > 0 ? M[i][j] : 0;
                        F[i][j] = F[i][j] > 0 ? F[i][j] : 0;
                        E[i][j] = E[i][j] > 0 ? E[i][j] : 0;

                        if( M[i][j] > maxScore ){ // please do not change > to >=, since we are doing local alignment
                                                        // >= will omit the first similar fragments
                                        maxScore = M[i][j];
                                        endPosition1=i; // this means the endPosition1 and endPosition2 is 1 based coordinate
                                        endPosition2=j;
                        }
                    }
            }
            ++a;
            --b;
         }

    }
}

int main()
{
    std::string dna1 = "AGCTGACCTGAACCTGAGAGCTGACCAGGCTGAAGGCTGACCTGGAAGCTGAGCTGAAGGCTGACCTGAAGCTGAAGGCTGAAGCTGACCTGGAAGGCTGACCTGGAAGCTGAGCTGACCTGAAGGCTGAAGGCTGACCTGGAAGCTGACCTGGAAGGCTGAAGGCTGAAGCTGACCTGAAGCTGAAGGCTGACCTGGAAGCTGAGCTGACCTGAAGGCTGACCTGAAGGCTGAAGCTGACCTGGAAGCTGAAGCTGAAGGCTGACCTGGAAGCTGACCTGAAGGCTGAAGGCTGACCTGGAAGCTGAGCTGACCTGAAGGCTGACCTGGAAGCTGACCTGAAGGCTGAAGCTGAAGGCTGAAGCTGACCTGGAAGGCTGACCTGGAAGCTGACCTGGAAGGCTGACCTGAAGCTGACCTGGAAGCTGAAGCTGAAGGCTGACCTGGAAGCTGAAGCTGAAGGCTGACCTGAAGCTGAAGGCTGACCTGGAAGCTGACCTGAAGCTGACCTGGAAGGCTGAAGCTGAAGCTGAAGGCTGAAGCTGACCTGGAAGGCTGACCTGAAGGCTGAAGCTGAAGGCTGACCTGGAAGGCTGAAGCTGAAGGCTGACCTGGAAGCTGACCTGAAGGCTGAAGCTGAAGGCTGACCTGGAAGGCTGAAGCTGAAGGCTGAAGCTGACCTGGAAGGCTGACCTGAAGCTGAAGGCTGAAGCTGAAGGCTGACCTGGAAGGCTGAAGCTGAAGGCTGACCTGGAAGCTGACCTGAAGGCTGAAGCTGAAGGCTGACCTGGAAGGCTGAAGCTGAAGGCTGACCTGGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGACCTGAAGCTGAAGGCTGACCTGGAAGGCTGAAGCTGAAGGCTGACCTGAAGCTGAAGGCTGAAGCTGACCTGGAAGGCTGAAGCTGAAGGCTGACCTGGAAGGCTGAAGCTGAAGGCTGACCTGGAAGGCTGAAGCTGAAGGCTGAAGCTGACCTGGAAGGCTGAAGCTGAAGGCTGACCTGGAAGGCTGAAGCTGAAGGCTGACCTGGAAGCTGACCTGAAGGCTGAAGCTGAAGGCTGACCTGGAAGGCTGAAGCTGAAGGCTGACCTGGAAGGCTGAAGCTGAAGGCTGACCTGGAAGGCTGAAGCTGAAGGCTGAAGCTGACCTGGAAGGCTGAAGCTGAAGGCTGACCTGGAAGGCTGAAGCTGAAGGCTGACCTGGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGA";
    std::string dna2 = "GCTGACC";

    int32_t length1 = dna1.length();
    int32_t length2 = dna2.length();

    int8_t *seq1 = new int8_t[length1];
    int8_t *seq2 = new int8_t[length2];

    for (int i = 0; i < length1; ++i)
        seq1[i] = dna1[i];
    for (int i = 0; i < length2; ++i)
        seq2[i] = dna2[i];

    int match_score = 2;
    int mismatch_score = -1;
    int open_gap_penalty = -2;
    int extend_gap_penalty = -1;

    int32_t maxScore, endPosition1, endPosition2;

    SmithWaterman(seq1, seq2, length1, length2, match_score, mismatch_score, open_gap_penalty, extend_gap_penalty, maxScore, endPosition1, endPosition2);

    std::cout << "Smith-Waterman Result:\n";
    std::cout << "Max Score: " << maxScore << "\n";
    std::cout << "End Position in Seq1: " << endPosition1 << "\n";
    std::cout << "End Position in Seq2: " << endPosition2 << "\n";

    SmithWaterman_Diag(seq1, seq2, length1, length2, match_score, mismatch_score, open_gap_penalty, extend_gap_penalty, maxScore, endPosition1, endPosition2);

    std::cout << "Smith-Waterman Diagonal Result:\n";
    std::cout << "Max Score: " << maxScore << "\n";
    std::cout << "End Position in Seq1: " << endPosition1 << "\n";
    std::cout << "End Position in Seq2: " << endPosition2 << "\n";

    delete[] seq1;
    delete[] seq2;

    return 0;
}