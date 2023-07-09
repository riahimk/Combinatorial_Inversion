# Combinatorial_Inversion
Combinatorial inversion is a new method that uses combinatorics to invert non-singular square matrices

Matrix Inversion Algorithms
Welcome to our repository! Here, we host novel and efficient matrix inversion algorithms developed as part of our research. These algorithms emphasize on triangular matrices and leverage combinatorial techniques and Strassen's fast matrix multiplication method. Our approach allows for a fully parallelizable execution, paving the way for efficient implementations on parallel computing architectures.

Featured Algorithms
CRIT: Column Recursive inverse of triangular matrices.
COMBRIT: Combinatorial (Block) Recursive inverse of triangular matrices.
SKUL: Inverse factorization with augmented classical $LU$ factorization.
SQR: Inverse factorization with augmented classical $QR$ factorization.
BRSI: Inverse Factorization based on recursive, split, and block inverse of triangular matrices.
Key Features
Efficient Matrix Inversions: Our algorithms offer highly efficient and direct computation of the inverse of triangular matrices.
Combines Multiple Techniques: We utilize combinatorial calculations, matrix splitting, and recurrent formalism for improved performance.
Enhanced Preconditioning: The algorithms have immense potential as preconditioners to accelerate Krylov subspace iterative methods and address large-scale systems of linear equations more efficiently.
Validation and Results
The proposed algorithms have been rigorously tested and validated. Our numerical tests demonstrate that these algorithms outperform traditional techniques, especially with larger matrices.

How to Use
You can download or clone this repository to use the algorithms in your own code. Each algorithm is implemented in a separate MATLAB script file.

bash
Copy code
git clone https://github.com/yourgithubusername/yourrepositoryname.git
For usage instructions, please refer to the individual documentation in each script file.

Further Research
While our research and implementations demonstrate the practicality of our proposed algorithms, we continue to explore advanced matrix inversion techniques. We aim to pave the way for improved numerical linear algebra algorithms and the development of effective preconditioners for various applications.


