all:
	mex ~/thesis/cauchy_lib/box_corr_lib.c -c
	mex ~/thesis/cauchy_lib/cauchy_swarz_ncc_match.c -c
	mex   mex_cauchy_ncc.c box_corr_lib.o  cauchy_swarz_ncc_match.o 
	mex   mex_cauchy_ncc_linear_idx.c box_corr_lib.o  cauchy_swarz_ncc_match.o 
	mex   mex_cauchy_ncc_linear_idx_clean.c box_corr_lib.o  cauchy_swarz_ncc_match.o 
	rm *.o -rf
