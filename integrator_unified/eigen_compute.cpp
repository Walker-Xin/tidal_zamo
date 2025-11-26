/**
 * @brief Computes eigenvalues and eigenvectors of a tidal tensor from trajectory data
 *
 * This function reads trajectory data from an input file, computes the tidal tensor
 * and its eigenvalues/eigenvectors at each point, and writes the results to an output file.
 *
 * @param filename_i Input filename containing trajectory data with format:
 *                   t r chi phi dt dr dchi dphi tau
 *
 * @note Output file is created in "data/" directory with name format:
 *       "eigensystem_a_[spin]_Q_[charge]_lambda_[lam]_[pro/ret].dat"
 */
void eigen_compute(char filename_i[128])
{
    // Variable declarations
    double t, r, chi, theta, phi;           // Complex variables for position and angles
    double dt, dr, dchi, dtheta, dphi, tau; // Complex variables for derivatives

    complex<double> C[4][4];        // 4x4 complex matrix
    Eigen::Matrix4cf C_eigen(4, 4); // Eigen matrix for eigenvalue computation

    int linenumber; // Line number counter

    char filename_o[128]; // Output filename
    FILE *foutput;        // Output file pointer

    // Open data file
    fstream newfile;
    newfile.open(filename_i, ios::in); // Open input file in read mode

    // Open output file
    if (retrograde == 0)
    {
        sprintf(filename_o, "data/eigensystem_a_%.2f_Q_%.2f_lambda_%.2f_pro.dat", spin, charge, lam); // Format output filename
    }
    else
    {
        sprintf(filename_o, "data/eigensystem_a_%.2f_Q_%.2f_lambda_%.2f_ret.dat", spin, charge, lam); // Format output filename
    }

    foutput = fopen(filename_o, "w"); // Open output file in write mode

    linenumber = 0; // Initialize line number counter

    if (newfile.is_open())
    {
        string tp; // Temporary string to hold each line of the file
        while (getline(newfile, tp))
        {
            // Extract data from tp, separated by space
            istringstream data(tp);
            data >> t >> r >> chi >> phi >> dt >> dr >> dchi >> dphi >> tau;

            // Calculate theta and dtheta from chi = cos(theta)
            theta = acos(chi);
            dtheta = -dchi / sqrt(1 - chi * chi);

            // Compute tidal tensor
            if (spherical == 0)
            {
                ClocH(r, theta, dt, dr, dtheta, dphi, C);
            }
            else
            {
                ClocH_spherical(r, theta, dt, dr, dtheta, dphi, C);
            }

            // Calculate eigenvalues and eigenvectors
            C_eigen = conversion(C);                                  // Convert C to an Eigen matrix
            Eigen::ComplexEigenSolver<Eigen::Matrix4cf> ces(C_eigen); // initialise eigen solver

            // Write four eigenvalues in one line in the output file
            fprintf(foutput, "%e %e %e %e ", ces.eigenvalues().col(0)[0].real(), ces.eigenvalues().col(0)[1].real(), ces.eigenvalues().col(0)[2].real(), ces.eigenvalues().col(0)[3].real());

            // Write four eigenvectors in four size-4 tuples in one line in the output file
            fprintf(foutput, "(%e,%e,%e,%e) (%e,%e,%e,%e) (%e,%e,%e,%e) (%e,%e,%e,%e)\n",
                    ces.eigenvectors().col(0)[0].real(), ces.eigenvectors().col(0)[1].real(), ces.eigenvectors().col(0)[2].real(), ces.eigenvectors().col(0)[3].real(),
                    ces.eigenvectors().col(1)[0].real(), ces.eigenvectors().col(1)[1].real(), ces.eigenvectors().col(1)[2].real(), ces.eigenvectors().col(1)[3].real(),
                    ces.eigenvectors().col(2)[0].real(), ces.eigenvectors().col(2)[1].real(), ces.eigenvectors().col(2)[2].real(), ces.eigenvectors().col(2)[3].real(),
                    ces.eigenvectors().col(3)[0].real(), ces.eigenvectors().col(3)[1].real(), ces.eigenvectors().col(3)[2].real(), ces.eigenvectors().col(3)[3].real());

            linenumber++;
        }
        newfile.close(); // Close input file
    }

    // Close output file
    fclose(foutput);
    printf("Number of lines: %d\n", linenumber); // Print number of lines processed
}