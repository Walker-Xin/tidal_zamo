void metric(double r, double chi, double g[4][4])
{
    double t2 = pow(r, 2);
    double t3 = pow(chi, 2);
    double t5 = pow(spin, 2);
    double t6 = t3 * t5;
    double t9 = pow(charge, 2);
    double t10 = -2 * M * r;
    double t7 = t2 + t6;
    double t8 = 1 / t7;
    double t13 = -1 + t3;
    double t14 = t10 + t9;
    double t15 = -(spin * t13 * t14 * t8);
    double t19 = -t3;
    double t20 = 1 + t19;
    double t16 = t10 + t2 + t5 + t9;
    g[0][0] = -(t8 * (t10 + t2 + t6 + t9));
    g[0][1] = 0;
    g[0][2] = 0;
    g[0][3] = t15;
    g[1][0] = 0;
    g[1][1] = t7 / t16;
    g[1][2] = 0;
    g[1][3] = 0;
    g[2][0] = 0;
    g[2][1] = 0;
    g[2][2] = t7 / t20;
    g[2][3] = 0;
    g[3][0] = t15;
    g[3][1] = 0;
    g[3][2] = 0;
    g[3][3] = t20 * (t13 * t16 * t5 + pow(t2 + t5, 2)) * t8;
}

void metric_spherical(double r, double chi, double g[4][4])
{
    double t6 = pow(r, 2);
    double t3 = pow(charge, 2);
    double t5 = -2 * M * r;
    double t7 = t3 + t5 + t6;
    double t11 = pow(chi, 2);
    g[0][0] = -(t7 / pow(r, 2));
    g[0][1] = 0;
    g[0][2] = 0;
    g[0][3] = 0;
    g[1][0] = 0;
    g[1][1] = t6 / t7;
    g[1][2] = 0;
    g[1][3] = 0;
    g[2][0] = 0;
    g[2][1] = 0;
    g[2][2] = t6 / (1 - t11);
    g[2][3] = 0;
    g[3][0] = 0;
    g[3][1] = 0;
    g[3][2] = 0;
    g[3][3] = -((-1 + t11) * t6);
}