#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <iostream>
#include <sys/time.h>
#include <math.h>
#include <omp.h>

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
                        int32_t &endPosition2) {

    int32_t **M = new int32_t *[length1 + 1];
    int32_t **E = new int32_t *[length1 + 1];
    int32_t **F = new int32_t *[length1 + 1];

    int32_t i, j;
    for (i = 0; i < (length1 + 1); ++i) {
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

    int32_t intervals = 200;
    int32_t rows_chunks = (length1 + intervals - 1) / intervals;
    int32_t cols_chunks = (length2 + intervals - 1) / intervals;

    int32_t local_maxScore = 0;
    int32_t local_endPosition1 = 0;
    int32_t local_endPosition2 = 0;
         for (int32_t diag = 0; diag < rows_chunks + cols_chunks - 1; ++diag) {
            #pragma omp parallel for num_threads(10) private(i, j) shared(M, E, F, seq1, seq2, maxScore, endPosition1, endPosition2)
            for(int32_t a = std::max(0, diag - cols_chunks + 1); a <= std::min(diag, rows_chunks) ; ++a){
                int32_t b = diag - a;
                for (i = a * intervals; i <= std::min((a + 1) * intervals, length1); ++i) {
                    for (j = b * intervals; j <= std::min((b + 1) * intervals, length2); ++j) {
                        if (i == 0 || j == 0) continue; // Skip processing for i=0 or j=0
                        E[i][j] = (_open_gap_penalty + M[i][j - 1]) > (_extend_gap_penalty + E[i][j - 1])
                                  ? (_open_gap_penalty + M[i][j - 1])
                                  : (_extend_gap_penalty + E[i][j - 1]);
                        F[i][j] = (_open_gap_penalty + M[i - 1][j]) > (_extend_gap_penalty + F[i - 1][j])
                                  ? (_open_gap_penalty + M[i - 1][j])
                                  : (_extend_gap_penalty + F[i - 1][j]);

                        int32_t match = (seq1[i - 1] == seq2[j - 1]) ? (_match_score) : (_mismatch_score);
                        match = match + M[i - 1][j - 1];

                        M[i][j] = (match > E[i][j]) ? (match) : E[i][j];
                        M[i][j] = M[i][j] > F[i][j] ? M[i][j] : F[i][j];
                        M[i][j] = M[i][j] > 0 ? M[i][j] : 0;
                        F[i][j] = F[i][j] > 0 ? F[i][j] : 0;
                        E[i][j] = E[i][j] > 0 ? E[i][j] : 0;

                        if (M[i][j] > local_maxScore) {
                            local_maxScore = M[i][j];
                            local_endPosition1 = i;
                            local_endPosition2 = j;
                        }
                    }
                }
            }
                if (local_maxScore > maxScore) {
                    maxScore = local_maxScore;
                    endPosition1 = local_endPosition1;
                    endPosition2 = local_endPosition2;
                }
    }

    // Memory cleanup
    for (i = 0; i <= length1; ++i) {
        delete[] M[i];
        delete[] E[i];
        delete[] F[i];
    }
    delete[] M;
    delete[] E;
    delete[] F;
}
int main()
{
    std::string dna1 = "AGCTGACCTGAACCTGAGAGCTGACCAGGCTGAAGGCTGACCTGGAAGCTGAGCTGAAGGCTGACCTGAAGCTGAAGGCTGAAGCTGACCTGGAAGGCTGACCTGGAAGCTGAGCTGACCTGAAGGCTGAAGGCTGACCTGGAAGCTGACCTGGAAGGCTGAAGGCTGAAGCTGACCTGAAGCTGAAGGCTGACCTGGAAGCTGAGCTGACCTGAAGGCTGACCTGAAGGCTGAAGCTGACCTGGAAGCTGAAGCTGAAGGCTGACCTGGAAGCTGACCTGAAGGCTGAAGGCTGACCTGGAAGCTGAGCTGACCTGAAGGCTGACCTGGAAGCTGACCTGAAGGCTGAAGCTGAAGGCTGAAGCTGACCTGGAAGGCTGACCTGGAAGCTGACCTGGAAGGCTGACCTGAAGCTGACCTGGAAGCTGAAGCTGAAGGCTGACCTGGAAGCTGAAGCTGAAGGCTGACCTGAAGCTGAAGGCTGACCTGGAAGCTGACCTGAAGCTGACCTGGAAGGCTGAAGCTGAAGCTGAAGGCTGAAGCTGACCTGGAAGGCTGACCTGAAGGCTGAAGCTGAAGGCTGACCTGGAAGGCTGAAGCTGAAGGCTGACCTGGAAGCTGACCTGAAGGCTGAAGCTGAAGGCTGACCTGGAAGGCTGAAGCTGAAGGCTGAAGCTGACCTGGAAGGCTGACCTGAAGCTGAAGGCTGAAGCTGAAGGCTGACCTGGAAGGCTGAAGCTGAAGGCTGACCTGGAAGCTGACCTGAAGGCTGAAGCTGAAGGCTGACCTGGAAGGCTGAAGCTGAAGGCTGACCTGGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGACCTGAAGCTGAAGGCTGACCTGGAAGGCTGAAGCTGAAGGCTGACCTGAAGCTGAAGGCTGAAGCTGACCTGGAAGGCTGAAGCTGAAGGCTGACCTGGAAGGCTGAAGCTGAAGGCTGACCTGGAAGGCTGAAGCTGAAGGCTGAAGCTGACCTGGAAGGCTGAAGCTGAAGGCTGACCTGGAAGGCTGAAGCTGAAGGCTGACCTGGAAGCTGACCTGAAGGCTGAAGCTGAAGGCTGACCTGGAAGGCTGAAGCTGAAGGCTGACCTGGAAGGCTGAAGCTGAAGGCTGACCTGGAAGGCTGAAGCTGAAGGCTGAAGCTGACCTGGAAGGCTGAAGCTGAAGGCTGACCTGGAAGGCTGAAGCTGAAGGCTGACCTGGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGA";
    std::string dna2 = "AGCTGACCTGAGCTGAAGCTGAAGCTGAAGCTGACCTGAAGCTGAAGCTGACCTGAAGCTGACCTGAGCTGAAGGCTGACCTGGAAGCTGACCTGAGCTGAAGGCTGAAGCTGACCTGAGCTGACCTGAAGGCTGAAGCTGAGCTGAAGGCTGACCTGGAAGGCTGAAGCTGAGCTGACCTGAAGGCTGACCTGAAGGCTGAAGCTGACCTGAAGGCTGAAGCTGACCTGGAAGCTGACCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGACCTGGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGACCTGGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGACCTGGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGACCTGGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGACCTGGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGACCTGGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGACCTGGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGAAGGCTGAAGCTGA";

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
    int32_t maxScore_diag, endPosition1_diag, endPosition2_diag;
    int w = 10;
    int h = 10;
    int32_t *maxScore_array = new int32_t[h];
    int32_t *endPosition1_array = new int32_t[h];
    int32_t *endPosition2_array = new int32_t[h];
    int32_t *maxScore_diag_array = new int32_t[h];
    int32_t *endPosition1_diag_array = new int32_t[h];
    int32_t *endPosition2_diag_array = new int32_t[h];

    long long sum0 = 0, sum1 = 0;
    struct timeval t0, t1;
    for (int r = 0; r < (w + h); ++r){
       gettimeofday(&t0, NULL);
       SmithWaterman(seq1, seq2, length1, length2, match_score, mismatch_score, open_gap_penalty, extend_gap_penalty, maxScore, endPosition1, endPosition2);
       gettimeofday(&t1, NULL);

       if(r >= w){sum0 += (t1.tv_sec - t0.tv_sec) * 1000000 + (t1.tv_usec - t0.tv_usec);
                  maxScore_array[r-w] = maxScore;
                  endPosition1_array[r-w] = endPosition1;
                  endPosition2_array[r-w] = endPosition2;
                  }
    }

    for (int r = 0; r < (w + h); ++r){
       gettimeofday(&t0, NULL);
       SmithWaterman_Diag(seq1, seq2, length1, length2, match_score, mismatch_score, open_gap_penalty, extend_gap_penalty, maxScore_diag, endPosition1_diag, endPosition2_diag);
       gettimeofday(&t1, NULL);

       if(r >= w){sum1 += (t1.tv_sec - t0.tv_sec) * 1000000 + (t1.tv_usec - t0.tv_usec);
                  maxScore_diag_array[r-w] = maxScore_diag;
                  endPosition1_diag_array[r-w] = endPosition1_diag;
                  endPosition2_diag_array[r-w] = endPosition2_diag;
                  }
    }

    int correct = 1;
    for (int i = 0; i < h; ++i){
        if (fabs(maxScore_array[i] - maxScore_diag_array[i]) >1e-6 ||
            fabs(endPosition1_array[i] - endPosition1_diag_array[i]) >1e-6 ||
            fabs(endPosition2_array[i] - endPosition2_diag_array[i]) >1e-6){
                 correct = 0;
                printf("%d\t%d\t%d\n", maxScore_array[i], endPosition1_array[i], endPosition2_array[i]);
                printf("%d\t%d\t%d\n", maxScore_diag_array[i], endPosition1_diag_array[i], endPosition2_diag_array[i]);
          }
    }
    printf("%d\t%lf\t%lf\n", correct, sum0/((double)(h*1.0)),sum1/((double)(h*1.0)));

    delete[] seq1;
    delete[] seq2;
    delete[] maxScore_array;
    delete[] endPosition1_array;
    delete[] endPosition2_array;
    delete[] maxScore_diag_array;
    delete[] endPosition1_diag_array;
    delete[] endPosition2_diag_array;

    return 0;
}
