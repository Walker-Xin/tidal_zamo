// Function to convert a 4x4 matrix of std::complex<double> to an Eigen::Matrix4cf
// Parameters:
// - C: A 4x4 matrix of std::complex<double>
// Returns:
// - An Eigen::Matrix4cf containing the same values as the input matrix
Eigen::Matrix4cf conversion(std::complex<double> C[4][4]) 
{
    // Convert the C matrix to a Matrix4cf
    Eigen::Matrix4cf C_eigen(4, 4);

    // Loop through each element of the 4x4 matrix
    for (int i = 0; i < 4; i++) 
    {
        for (int j = 0; j < 4; j++) 
        {
            // Assign the value from C to the corresponding position in C_eigen
            C_eigen(i, j) = C[i][j];
        }
    }

    // Return the converted Eigen matrix
    return C_eigen;
}