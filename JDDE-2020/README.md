Codes to accompany "Eigenvalues and delay differential equations: periodic coefficients, impulses and rigorous numerics".
Current version: 1.0

Contents: 
- folder FP : Floating-point implementation of -just- the monodromy operator discretization. If you only want to approximate eigenvalues and check stabiity of an impulsive DDE, this is what you want.
- zip file monodromy_2020.zip : convenient download folder containing all necessary MATLAB function files needed to complete computer-assisted proofs from the paper or monodromy operator discretization.

Current development:
- optimize by using multidimensional arrays for inputs (ans subsequent processing) rather than cell arrays
- vectorize some loops
