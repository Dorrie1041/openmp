#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <iostream>
#include <sys/time.h>
#include <math.h>
#include <omp.h>
#include <fstream>
#include <sstream>

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
                        int32_t &endPosition2,
                        int32_t &intervals,
                        int num_threads) {
    int32_t i, j;

    int32_t *E = new int32_t[length1 + 1];
    std::fill_n(E, length1 + 1, 0);

    int32_t *F = new int32_t[length2 + 1];
    std::fill_n(F, length2 + 1, 0);

    int32_t *C = new int32_t[length1 + 1];
    std::fill_n(C, length1 + 1, 0);

    int32_t **M = new int32_t *[2];
    for (i = 0; i < 2; ++i) {
        M[i] = new int32_t[length2 + 1];
        std::fill_n(M[i], length2 + 1, 0);
    }

    endPosition1 = 0;
    endPosition2 = 0;
    maxScore = 0;

    int32_t rows_chunks = (length1 + intervals - 1) / intervals;
    int32_t cols_chunks = (length2 + intervals - 1) / intervals;
    int32_t max_num = std::max(rows_chunks, cols_chunks);
    int32_t **R = new int32_t *[2];
    for (i = 0; i < 2; ++i) {
        R[i] = new int32_t[max_num + 1];
	std::fill_n(R[i], max_num + 1, 0);
    }
    for (int32_t diag = 0; diag < rows_chunks + cols_chunks - 1; ++diag) {
        int current_chunks = diag % 2;
        int previous_chunks = 1 - current_chunks;

	int32_t *temp_R = new int32_t[max_num + 1];
	std::fill_n(temp_R, max_num + 1, 0);
        #pragma omp parallel num_threads(num_threads) private(i, j)
        {
            int32_t localMaxScore = maxScore;
            int32_t localEndPosition1 = endPosition1;
            int32_t localEndPosition2 = endPosition2;

            #pragma omp for
            for (int32_t a = std::max(0, diag - cols_chunks + 1); a <= std::min(diag, rows_chunks); ++a) {
                int32_t b = diag - a;
                int32_t *temp_C = new int32_t[length1 + 1];
                std::copy(C, C + length1 + 1, temp_C);

                for (int32_t ia = a * intervals; ia < std::min((a + 1) * intervals, length1); ++ia) {
                    for (int32_t jb = b * intervals; jb < std::min((b + 1) * intervals, length2); ++jb) {
                        int i = ia + 1;
                        int j = jb + 1;
                        int current_i = i % 2;
                        int previous_i = 1 - current_i;

                        if (j == b * intervals + 1) {
                            E[i] = (_open_gap_penalty + temp_C[i]) > (_extend_gap_penalty + E[i])
                                       ? (_open_gap_penalty + temp_C[i])
                                       : (_extend_gap_penalty + E[i]);
                        } else {
                            E[i] = (_open_gap_penalty + M[current_i][j - 1]) > (_extend_gap_penalty + E[i])
                                       ? (_open_gap_penalty + M[current_i][j - 1])
                                       : (_extend_gap_penalty + E[i]);
                        }

                        F[j] = (_open_gap_penalty + M[previous_i][j]) > (_extend_gap_penalty + F[j])
                                   ? (_open_gap_penalty + M[previous_i][j])
                                   : (_extend_gap_penalty + F[j]);

                        int32_t match = (seq1[i - 1] == seq2[j - 1]) ? _match_score : _mismatch_score;
                        if (j == b * intervals + 1 && j != 1) {
                            if (i == a * intervals + 1 && i != 1) {
			        match += R[current_chunks][a - 1];
                            } else {
                                match += temp_C[i - 1];
                            }
                        } else {
                            match += M[previous_i][j - 1];
                        }

                        M[current_i][j] = (match > E[i]) ? (match) : E[i];
                        M[current_i][j] = M[current_i][j] > F[j] ? M[current_i][j] : F[j];
                        M[current_i][j] = M[current_i][j] > 0 ? M[current_i][j] : 0;
                        F[j] = F[j] > 0 ? F[j] : 0;
                        E[i] = E[i] > 0 ? E[i] : 0;
			
                        if (M[current_i][j] > localMaxScore) {
                            localMaxScore = M[current_i][j];
                            localEndPosition1 = i;
                            localEndPosition2 = j;
                        }

                        if (j == std::min((b + 1) * intervals, length2)) {
                            if (i == std::min((a + 1) * intervals, length1)) {
                               temp_R[a] = M[current_i][j];
                            }
			    
			    C[i] = M[current_i][j];
		       }
                    }
                }

                delete[] temp_C;
            }

            #pragma omp critical
            {
                if (localMaxScore > maxScore) {
                    maxScore = localMaxScore;
                    endPosition1 = localEndPosition1;
                    endPosition2 = localEndPosition2;
                }
		
		for (int i = 0; i < (max_num + 1); ++i){
		  R[current_chunks][i] = temp_R[i];
		}
            }
        }
	delete[] temp_R;
    }

    for (i = 0; i < 2; ++i) {
        delete[] M[i];
    }
    for (i = 0; i < 2; ++i) {
        delete[] R[i];
    }
    delete[] R;
    delete[] M;
    delete[] E;
    delete[] F;
    delete[] C;
}

int main(int argc, char **argv){
    if (argc < 4) {
        std::cerr << "Usage: " << argv[0] << " <interval> <threads> <file_path> " << std::endl;
        return 1;
    }
       int32_t interval = atoi(argv[1]);
       int threads = atoi(argv[2]);
       std::string file_path = argv[3];
       std::string dna1, dna2;
       std::ifstream file(file_path);
       if (!file.is_open()){
          std::cerr << "Error opening file: " << file_path << std::endl;
       }

       std::string line;
       bool found_reference = false;
       while (std::getline(file, line)){
           if (line.rfind(">reference", 0) == 0){
                  found_reference = true;
             } else if (line.rfind(">query", 0) == 0){
                  found_reference = false;
             } else if (found_reference){
                  dna1 += line;
             } else {
                  dna2 += line;
             }
        }

       file.close();
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
       int w = 5;
       int h = 5;
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
          SmithWaterman_Diag(seq1, seq2, length1, length2, match_score, mismatch_score, open_gap_penalty, extend_gap_penalty, maxScore_diag, endPosition1_diag, endPosition2_diag, interval, threads);
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
       printf("%d\t%d\t%d\t%d\t%d\t%lf\t%lf\n", length1, length2, interval, threads, correct, sum0/((double)(h*1.0)),sum1/((double)(h*1.0)));

       delete[] seq1;
       delete[] seq2;
       delete[] maxScore_array;
       delete[] endPosition1_array;
       delete[] endPosition2_array;
       delete[] maxScore_diag_array;
       delete[] endPosition1_diag_array;
       delete[] endPosition2_diag_array;

       }
