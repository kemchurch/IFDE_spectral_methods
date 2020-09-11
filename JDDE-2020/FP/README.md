Contents:
- FP.zip : All .m files, demo scripts and a manual file needed to run the primary implementation.
- M_truncation.m : Primary implementation of the Chebyshev spectral method that is applicable for arbitrary positive integer period and delay. Calls every other function in this directory EXCEPT E_truncation_pq1.m and M_truncation_pq1.m.
- M_truncation_pq1.m : Old implementation that applies only to the special case p=q=1. May be useful for educational purposes, but now defunct.