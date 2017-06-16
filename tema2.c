#include "utils.h"

#define MAX(a,b) (((a)>(b))?(a):(b))

/* Multiplication function
 */
double *f_dgemm(char transa, char transb,
	int M, int N, int K,
	double alpha, double *A, int lda,
	double *B, int ldb,
	double beta, double *C, int ldc) {

	/* macros used for matrix access in row-major order */
	#define A(i, j) A[i + j * lda]
	#define B(i, j) B[i + j * ldb]
	#define C(i, j) C[i + j * ldc]

	int i, j, l;
	int no_trans_a = 0, no_trans_b = 0;
	int nrowa, nrowb;

	no_trans_a = (transa == 'N');
	no_trans_b = (transb == 'N');

	/* sanity checks */
	nrowa = K;

    if (no_trans_a == 1) {
    	nrowa = M;
    }

    nrowb = N;

    if (no_trans_b == 1)
    	nrowb = K;

    if (no_trans_a == 0 && transa != 'C' && transa != 'T')
    	error("Wrong transa.");
    else if (no_trans_b == 0 && transb != 'C' && transb != 'T')
    	error("Wrong transb.");
    else if (M < 0)
    	error("M < 0");
    else if (N < 0)
    	error("N < 0");
    else if (K < 0)
    	error("K < 0");
    else if (lda < MAX(1, nrowa))
    	error("lda < MAX(1, nrowa).");
    else if (ldb < MAX(1, nrowb))
    	error("ldb < MAX(1, nrowb).");
    else if (ldc < MAX(1, M))
    	error("ldc < MAX(1, M).");


    /* quick return */
	if ((M == 0 || N == 0 || (double_equal(alpha, 0.0) ||
		(K == 0))) && (double_equal(beta, 1.0)))
		return C;

	/* I use pointers for memory access */
	/* do less operations when add / multiply: 2 adds 
	 * and 1 mul instead of 2N adds and N multiplies
	 */
	/* I traverse the matrices in row-major order
	 * because I read the transpose directly in main
	 */
	if (double_equal(alpha, 0.0)) {
		if (double_equal(beta, 0.0)) {
			for (j = 0; j < N; ++j) {
				double *pc = &(C(0, j));
				for (i = 0; i < M; ++i) {
					*pc = 0.0;
					pc++;
				}
			}
		} else {
			for (j = 0; j < N; ++j) {
				double *pc = &(C(0, j));
				for (i = 0; i < M; ++i) {
					*pc = beta * (*pc);
					pc++;
				}
			}
		}
		return C;
	}

	/* also, I use a registers for the constants in the loops */
	/* all accesses are sequential, so we don't have to access 
	 * elements at addresses that are multiple of "size"
	 */
	if (no_trans_b == 1) {
		if (no_trans_a == 1) {
			for (j = 0; j < N; ++j) {
				if (double_equal(beta, 0.0) == 1) {
					double *pc = &(C(0, j));
					for (i = 0; i < M; ++i) {
						*pc = 0.0;
						pc++;
					}
				} else if (double_equal(beta, 1.0) == 0) {
					double *pc = &(C(0, j));
					for (i = 0; i < M; ++i) {
						*pc = beta * (*pc);
						pc++;
					}
				}

				double *pb = &(B(0, j));
				for (l = 0; l < K; ++l) {
					register double temp = alpha * (*pb);
					pb++;

					double *pa = &(A(0, l));
					double *pc = &(C(0, j));
					for (i = 0; i < M; ++i) {
						*pc += temp * (*pa);
						pc++;
						pa++;
					}
				}
			}
		} else if (no_trans_a == 0) {
			for (j = 0; j < N; ++j) {
				double *pc = &(C(0, j));
				for (i = 0; i < M; ++i) {

					double *pb = &(B(0, j));
					double *pa = &(A(0, i));
					register double temp = 0.0;
					for (l = 0; l < K; ++l) {
						temp += (*pa) * (*pb);
						pa++;
						pb++;
					}
					if (double_equal(beta, 0.0) == 1) {
						*pc = alpha * temp;
						pc++;
					} else {
						*pc = alpha * temp + beta * (*pc);
						pc++;
					}
				}
			}
		}
	} else if (no_trans_b == 0) {
		if (no_trans_a == 1) {
			for (j = 0; j < N; ++j) {
				if (double_equal(beta, 0.0) == 1) {
					double *pc = &(C(0, j));
					for (i = 0; i < M; ++i) {
						*pc = 0.0;
						pc++;
					}
				} else if (double_equal(beta, 1.0) == 0) {
					double *pc = &(C(0, j));
					for (i = 0; i < M; ++i) {
						*pc = beta * (*pc);
						pc++;
					}
				}

				double *pb = &(B(0, j));
				for (l = 0; l < K; ++l) {
					register double temp = alpha * (*pb);
					pb++;

					double *pa = &(A(0, l));
					double *pc = &(C(0, j));
					for (i = 0; i < M; ++i) {
						(*pc) += temp * (*pa);
						pa++;
						pc++;
					}
				}
			}
		} else if (no_trans_a == 0) {
			for (j = 0; j < N; ++j) {
				double *pc = &(C(0, j));
				for (i = 0; i < M; ++i) {

					double *pb = &(B(0, j));
					double *pa = &(A(0, i));
					register double temp = 0.0;

					for (l = 0; l < K; ++l) {
						temp += (*pa) * (*pb);
						pb++;
						pa++;
					}
					if (double_equal(beta, 0.0) == 1) {
						*pc = alpha * temp;
						pc++;
					} else {
						*pc = alpha * temp + beta * (*pc);
						pc++;
					}
				}
			}
		}
	}

	/* return the result */
	return C;
}

double my_atof(const char *value, int length) {
	int digit, found_point = 0;
	double double_val = 0, factor = 1;
	int i = 0;

	while (i < length) {

		/* if we meet the point, then, at the end we will divide the result
		* by the number of digits after this char(.)
		*/
		if (*value == '.') {
			found_point = 1;
			i++;
			value++;
			continue;
		}

		/* get a digit [0-9]*/
		digit = *value - '0';

		if (digit >= 0 && digit <= 9) {

			/* divide factor by 10 for every digit after the "."" */
			if (found_point)
				factor /= 10.0;

			/* compute the number */
			double_val = double_val * 10.0 + (double) digit;
		}

		i++;
		value++;
	}

	/* divide the whole number to the factor to obtain the double */
	return double_val * factor;
}


int main(int argc, char **argv)
{
	struct test **tests;
	int i, k, m, n;
	char path[256];
	void *AA, *BB, *CC;
	struct stat s;
	int fd, rc;
	int status;
    size_t file_size;
    char *mapped;
    FILE *file = NULL;

	memset(path, 0, 256);

	tests = (struct test **)malloc(sizeof(struct test*));
	(*tests) = (struct test *)calloc(MAXTESTS, sizeof(struct test));

	parse_config("tema2.cfg", tests);

	/* create out directory if necessary */
	struct stat st = {0};

	if (stat("./out", &st) == -1) {
    	mkdir("./out", 0700);
	}

	for (i = 0; i < MAXTESTS; i++) {

		struct test *my_test = &(*tests)[i];

		if (my_test->name[0] == '\0')
		    break;

		/* alloc memory for matrices using minimum padding */
		AA = (double *) malloc(my_test->M * my_test->K * sizeof(double) + 31);

		my_test->A = (double *)(((uintptr_t)AA + 31) & ~31);

		BB = (double *) malloc(my_test->K * my_test->N * sizeof(double) + 31);

		my_test->B = (double *)(((uintptr_t)BB + 31) & ~31);

		CC = (double *) malloc(my_test->M * my_test->N * sizeof(double) + 31);

		my_test->C = (double *)(((uintptr_t)CC + 31) & ~31);

		/* I read the transposed matrices, because I need to access them in f_dgemm in row-major order */

		/* read matrix A from file using mmap */

		strcpy(path, "input/");
		strcat(path, my_test->name);
		strcat(path, "_A.in");

		fd = open(path, O_RDONLY);
		if (fd < 0)
			error("Opening the file.");

		/* get the struct info */
		status = fstat(fd, &s);
		if (status < 0)
			error("File status.");

		file_size = s.st_size;

		mapped = mmap(0, file_size, PROT_READ, MAP_PRIVATE, fd, 0);
		if (mapped == MAP_FAILED)
			error("Mmap failed");

		memset(path, 0, 256);

		char *ptr;

		/* use a temp pointer for reading */ 
		ptr = mapped;
		char value[16];
		memset(value, 0, 16);

		for (m = 0; m < my_test->M; ++m) {
			for (k = 0; k < my_test->K; ++k) {

				/* read a value from file in a char[] */
				int ii = 0;
				while(*ptr != ' ') {
					if (*ptr != '\n') {
						value[ii] = *ptr;
						ii++;
					}
					ptr++;
				}
				/* skip the space */
				ptr++;

				/* convert the string to double and save it in the matrix */
				my_test->A[m + k * my_test->M] = my_atof(value, ii);

				memset(value, 0, 16);
			}
		}

		/* free the mapping */
		rc = munmap((void *)mapped, file_size);
		if (rc < 0)
			error("Munmap failed.");

		close(fd);

		/* did the same thing for B and C as above, too */

		/* read matrix B from file using mmap */
		strcpy(path, "input/");
		strcat(path, my_test->name);
		strcat(path, "_B.in");

		fd = open(path, O_RDONLY);
		if (fd < 0)
			error("Opening the file.");

		status = fstat(fd, &s);
		if (status < 0)
			error("File status.");

		file_size = s.st_size;

		mapped = mmap(0, file_size, PROT_READ, MAP_PRIVATE, fd, 0);
		if (mapped == MAP_FAILED)
			error("Mmap failed");

		memset(path, 0, 256);

		ptr = mapped;
		memset(value, 0, 16);

		/* here I read B or B**T, just because i want to access it sequentialy in the f_dgemm function */
		if (my_test->transb == 'N') {
			for (k = 0; k < my_test->K; ++k) {
				for (n = 0; n < my_test->N; ++n) {

					int ii = 0;
					while(*ptr != ' ') {
						if (*ptr != '\n') {
							value[ii] = *ptr;
							ii++;
						}
						ptr++;
					}
					ptr++;

					my_test->B[k + n * my_test->K] = my_atof(value, ii);

					memset(value, 0, 16);
				}
			}
		} else if (my_test->transb == 'T') {
			for (k = 0; k < my_test->K; ++k) {
				for (n = 0; n < my_test->N; ++n) {

					int ii = 0;
					while(*ptr != ' ') {
						if (*ptr != '\n') {
							value[ii] = *ptr;
							ii++;
						}
						ptr++;
					}
					ptr++;

					my_test->B[n + k * my_test->N] = my_atof(value, ii);

					memset(value, 0, 16);
				}
			}
		}

		rc = munmap((void *)mapped, file_size);
		if (rc < 0)
			error("Munmap failed.");

		close(fd);

		/* read matrix C from file using mmap */

		strcpy(path, "input/");
		strcat(path, my_test->name);
		strcat(path, "_C.in");

		fd = open(path, O_RDONLY);
		if (fd < 0)
			error("Opening the file.");

		status = fstat(fd, &s);
		if (status < 0)
			error("File status.");

		file_size = s.st_size;

		mapped = mmap(0, file_size, PROT_READ, MAP_PRIVATE, fd, 0);
		if (mapped == MAP_FAILED)
			error("Mmap failed");

		memset(path, 0, 256);

		ptr = mapped;
		memset(value, 0, 16);
		for (m = 0; m < my_test->M; ++m) {
			for (n = 0; n < my_test->N; ++n) {

				int ii = 0;
				memset(value, 0, 16);
				while(*ptr != ' ') {
					if (*ptr != '\n') {
						value[ii] = *ptr;
						ii++;
					}
					ptr++;
				}
				ptr++;

				my_test->C[m + n * my_test->M] = my_atof(value, ii);
			}
		}

		rc = munmap((void *)mapped, file_size);
		if (rc < 0)
			error("Munmap failed.");

		close(fd);

		/* execute the algorithm */
		my_test->C = f_dgemm(my_test->transa, my_test->transb, my_test->M, my_test->N,
			my_test->K, my_test->alpha, my_test->A, my_test->lda,
			my_test->B, my_test->ldb, my_test->beta, my_test->C, my_test->ldc);


		/* write the output to the file, each value having exactly 3 decimals*/
		strcpy(path, "out/");
		strcat(path, my_test->name);
		strcat(path, ".out");

		file = fopen(path, "w");

		memset(path, 0, 256);

		for (m = 0; m < my_test->M; ++m) {
			for (n = 0; n < my_test->N; ++n) {
				if (!fprintf(file, "%.3lf ", my_test->C[m + n * my_test->M]))
					break;
			}
			if (!fprintf(file, "%s", "\n"))
					break;
		}

		fclose(file);


		/* free memory */
		free(AA);
		free(BB);
		free(CC);
	}

	return 0;
}
