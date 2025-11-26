/**
 * @brief Performs ray tracing simulation for photon trajectories in a curved spacetime
 *
 * This function implements a Dormand-Prince 8(5,3) integrator to trace photon paths
 * in either spherical or non-spherical coordinate systems. It integrates the geodesic
 * equations of motion while monitoring various conserved quantities.
 *
 * The integration continues until one of these conditions is met:
 * - Photon crosses the singularity (r < 0)
 * - Photon escapes (r > 1.5 * initial radius)
 * - Maximum iterations reached
 *
 * During integration, it checks:
 * - Conservation of invariant length
 * - Conservation of energy
 * - Conservation of angular momentum
 * - Crossing of IBSO (innermost bound spherical orbit)
 *
 * Results are written to an output file in the format:
 * t r chi phi kt kr kchi kphi affine_parameter
 *
 * @param filename_o Output file path where trajectory data will be written
 * @param total_iterations Maximum number of integration steps to perform
 * @param error_tolerance Desired accuracy of the integration (default step error tolerance)
 *
 * @note Uses global variables:
 *  - spherical: Coordinate system flag (0 for non-spherical, 1 for spherical)
 *  - psi: Initial angle parameter
 *  - lz: Angular momentum
 *  - energy: Energy
 *  - x: IBSO radius
 */
void raytrace(char filename_o[128], int total_iterations, double error_tolerance)
{
    double denom;
    double t, r, th, chi, phi;
    double t0, r0, chi0, th0, phi0;
    double kt, kr, kth, kchi, kphi;
    double kt0, kr0, kth0, kchi0, kphi0;
    double h, hnext, aff;

    double g[4][4]; // metric tensor
    double diffs[8], diffs_new[8], vars[8], dvars[8], vars_temp[8], vars_out[8], vars_4th[8], vars_5th[8], k1[8], k2[8], k3[8], k4[8], k5[8], k6[8], k7[8], k8[8], k9[8], k10[8];
    double err[8], err1[8];
    double errtol, errmin, errmax;
    double safety, max_next, min_next;
    double temp, errscale, errval, errval1;
    double test;

    int errcheck = 0, crosscheck = 0;
    int stop_integration = 0;
    int counter;
    int length_check = 0, energy_check = 0, angmom_check = 0;

    FILE *foutput;

    // Set error tolerances and relevant variables
    errtol = error_tolerance;
    errmin = errtol / 10; // minimum error tolerance
    errmax = errtol * 10; // maximum error tolerance
    safety = 0.95;        // safety factor
    min_next = 0.333;     // minimum scaling factor
    max_next = 6;         // maximum scaling factor

    h = 1.0; // initial step size

    // Set initial conditions
    // Initial coordinates
    t0 = 0.0;
    r0 = 1.0e+6;
    th0 = psi;
    chi0 = cos(psi);
    phi0 = 0; // set as zero since spacetime is axisymmetric

    // Compute initial metric tensor
    if (spherical == 0)
    {
        metric(r0, chi0, g);
    }
    else
    {
        metric_spherical(r0, chi0, g);
    }

    // Initial velocities
    kth0 = 0;
    kchi0 = 0;

    if (spherical == 0)
    {
        denom = g[3][3] * g[0][0] - g[0][3] * g[0][3]; // used for convenience
        kt0 = -(lz * g[0][3] + energy * g[3][3]) / denom;
        kphi0 = (lz * g[0][0] + energy * g[0][3]) / denom;

        // Get initial r velocity via invariant length
        kr0 = sqrt(-1 - g[0][0] * kt0 * kt0 - g[2][2] * kth0 * kth0 - g[3][3] * kphi0 * kphi0 - 2 * g[0][3] * kt0 * kphi0) / sqrt(g[1][1]);
        kr0 = -kr0; // negative sign for inward going particle
    }
    else
    {
        kt0 = -energy / g[0][0];
        kphi0 = lz / g[3][3];
        // Get initial r velocity via invariant length
        kr0 = sqrt(-1 - g[0][0] * kt0 * kt0 - g[2][2] * kth0 * kth0 - g[3][3] * kphi0 * kphi0) / sqrt(g[1][1]);
        kr0 = -kr0; // negative sign for inward going particle
    }

    // Open output file
    foutput = fopen(filename_o, "w");

    // === DORMAND-PRINCE 8(5,3) INTEGRATOR ===
    // Initialise variables
    t = t0;
    r = r0;
    th = th0;
    chi = chi0;
    phi = phi0;
    aff = 0.0;

    kt = kt0;
    kr = kr0;
    kth = kth0;
    kchi = kchi0;
    kphi = kphi0;

    counter = 0;

    do
    {
        vars[0] = t;
        vars[1] = r;
        vars[2] = chi;
        vars[3] = phi;
        vars[4] = kt;
        vars[5] = kr;
        vars[6] = kchi;
        vars[7] = kphi;

        do
        {
            errcheck = 0;

            // Evaluation 1
            diffeqs(vars, diffs);
            for (int i = 0; i <= 7; i++)
            {
                k1[i] = h * diffs[i];
                vars_temp[i] = vars[i] + a21 * k1[i];
            }

            // Evaluation 2
            diffeqs(vars_temp, diffs);
            for (int i = 0; i <= 7; i++)
            {
                k2[i] = h * diffs[i];
                vars_temp[i] = vars[i] + a31 * k1[i] + a32 * k2[i];
            }

            // Evaluation 3
            diffeqs(vars_temp, diffs);
            for (int i = 0; i <= 7; i++)
            {
                k3[i] = h * diffs[i];
                vars_temp[i] = vars[i] + a41 * k1[i] + a43 * k3[i];
            }

            // Evaluation 4
            diffeqs(vars_temp, diffs);
            for (int i = 0; i <= 7; i++)
            {
                k4[i] = h * diffs[i];
                vars_temp[i] = vars[i] + a51 * k1[i] + a53 * k3[i] + a54 * k4[i];
            }

            // Evaluation 5
            diffeqs(vars_temp, diffs);
            for (int i = 0; i <= 7; i++)
            {
                k5[i] = h * diffs[i];
                vars_temp[i] = vars[i] + a61 * k1[i] + a64 * k4[i] + a65 * k5[i];
            }

            // Evaluation 6
            diffeqs(vars_temp, diffs);
            for (int i = 0; i <= 7; i++)
            {
                k6[i] = h * diffs[i];
                vars_temp[i] = vars[i] + a71 * k1[i] + a74 * k4[i] + a75 * k5[i] + a76 * k6[i];
            }

            // Evaluation 7
            diffeqs(vars_temp, diffs);
            for (int i = 0; i <= 7; i++)
            {
                k7[i] = h * diffs[i];
                vars_temp[i] = vars[i] + a81 * k1[i] + a84 * k4[i] + a85 * k5[i] + a86 * k6[i] + a87 * k7[i];
            }

            // Evaluation 8
            diffeqs(vars_temp, diffs);
            for (int i = 0; i <= 7; i++)
            {
                k8[i] = h * diffs[i];
                vars_temp[i] = vars[i] + a91 * k1[i] + a94 * k4[i] + a95 * k5[i] + a96 * k6[i] + a97 * k7[i] + a98 * k8[i];
            }

            // Evaluation 9
            diffeqs(vars_temp, diffs);
            for (int i = 0; i <= 7; i++)
            {
                k9[i] = h * diffs[i];
                vars_temp[i] = vars[i] + a101 * k1[i] + a104 * k4[i] + a105 * k5[i] + a106 * k6[i] + a107 * k7[i] + a108 * k8[i] + a109 * k9[i];
            }

            // Evaluation 10
            diffeqs(vars_temp, diffs);
            for (int i = 0; i <= 7; i++)
            {
                k10[i] = h * diffs[i];
                vars_temp[i] = vars[i] + a111 * k1[i] + a114 * k4[i] + a115 * k5[i] + a116 * k6[i] + a117 * k7[i] + a118 * k8[i] + a119 * k9[i] + a1110 * k10[i];
            }

            // Evaluation 11
            diffeqs(vars_temp, diffs);
            for (int i = 0; i <= 7; i++)
            {
                k2[i] = h * diffs[i];
                vars_temp[i] = vars[i] + a121 * k1[i] + a124 * k4[i] + a125 * k5[i] + a126 * k6[i] + a127 * k7[i] + a128 * k8[i] + a129 * k9[i] + a1210 * k10[i] + a1211 * k2[i];
            }

            // Evaluation 12
            diffeqs(vars_temp, diffs);
            for (int i = 0; i <= 7; i++)
            {
                k3[i] = h * diffs[i];
                k4[i] = b1 * k1[i] + b6 * k6[i] + b7 * k7[i] + b8 * k8[i] + b9 * k9[i] + b10 * k10[i] + b11 * k2[i] + b12 * k3[i];
                vars_out[i] = vars[i] + k4[i];
            }

            // Get two error estimates
            for (int i = 0; i <= 7; i++)
            {
                // Divide by h compared to the original code
                err[i] = (k4[i] - bhh1 * k1[i] - bhh2 * k9[i] - bhh3 * k3[i]) / h;
                err1[i] = (er1 * k1[i] + er6 * k6[i] + er7 * k7[i] + er8 * k8[i] + er9 * k9[i] + er10 * k10[i] + er11 * k2[i] + er12 * k3[i]) / h;
            }

            // Compute error
            errval = 0.0;
            for (int i = 0; i <= 7; i++)
            {
                errscale = errtol + errtol * max(fabs(vars[i]), fabs(vars_out[i]));
                errval += err[i] * err[i] / (errscale * errscale);
                errval1 += err1[i] * err1[i] / (errscale * errscale);
            }

            temp = errval + 0.01 * errval1;
            if (temp <= 0)
                temp = 1.0;

            errval = fabs(h) * errval * sqrt(1.0 / (8.0 * temp));

            // Adjust step size
            if (errval > 1) // step error too large
            {
                hnext = safety * pow(1 / errval, 0.125) * h;

                if (hnext < min_next * h)
                {
                    hnext = min_next * h; // do not allow step to decrease too much
                }

                errcheck = 1;
            }
            else // step successful
            {
                hnext = safety * pow(1 / errval, 0.125) * h;

                if (hnext > max_next * h)
                {
                    hnext = max_next * h; // do not allow step to increase too much
                }
                if (hnext < min_next * h)
                {
                    hnext = min_next * h; // do not allow step to decrease too much
                }

                errcheck = 0;
            }

            h = hnext;

        } while (errcheck == 1);

        // Error satisfactory, update values
        for (int i = 0; i <= 7; i++)
            vars[i] = vars_out[i];

        t = vars[0];
        r = vars[1];
        chi = vars[2];
        phi = vars[3];
        kt = vars[4];
        kr = vars[5];
        kchi = vars[6];
        kphi = vars[7];
        aff += h;

        // Record position
        fprintf(foutput, "%.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f\n", double(t), double(r), double(chi), double(phi), double(kt), double(kr), double(kchi), double(kphi), double(aff));

        counter++;

        // Check invariant length
        if (spherical == 0)
        {
            metric(r, chi, g);
            denom = g[3][3] * g[0][0] - g[0][3] * g[0][3]; // used for convenience
            kt = -(lz * g[0][3] + energy * g[3][3]) / denom;
            kphi = (lz * g[0][0] + energy * g[0][3]) / denom;
        }
        else
        {
            metric_spherical(r, chi, g);
            kt = -energy / g[0][0];
            kphi = lz / g[3][3];
            ;
        }
        test = 1 + g[0][0] * kt * kt + g[1][1] * kr * kr + g[2][2] * kchi * kchi + g[3][3] * kphi * kphi + 2 * g[0][3] * kt * kphi;

        if (fabs(test) > 1.0e-7 && length_check == 0)
        {
            printf("Invariant length not satisfied, iteration = %d\n, test = %f\n", counter, test);
            length_check = 1;
            // break;
        }

        // Check energy
        if (fabs(energy + g[0][0] * kt + g[0][3] * kphi) > 1.0e-7 && energy_check == 0)
        {
            printf("Energy not conserved, iteration = %d\n", counter);
            energy_check = 1;
            // break;
        }

        // Check angular momentum
        if (fabs(lz - g[3][3] * kphi - g[0][3] * kt) > 1.0e-7 && angmom_check == 0)
        {
            printf("Angular momentum not conserved, iteration = %d\n", counter);
            angmom_check = 1;
            // break;
        }

        // Check if r is less than x
        if (vars[1] < x && crosscheck == 0)
        {
            printf("Photon crossed IBSO at r = %f, iteration = %d\n", vars[1], counter);
            crosscheck = 1;
        }

        // Check if r is less than 0
        if (vars[1] < 0)
        {
            printf("Photon crossed singularity at r = %f, iteration = %d\n", vars[1], counter);
            break;
        }

        // Check if particle has escaped
        if (vars[1] > (1.5 * r0))
        {
            printf("Photon escaped at r = %f, iteration = %d\n", vars[1], counter);
            break;
        }

    } while (counter < total_iterations);
}
