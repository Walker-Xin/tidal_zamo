void diffeqs(double vars[], double diffs[])
{
	double t, r, chi, phi;
	double dt, dr, dchi, dphi;
	double christ[4][4][4];

	t = vars[0];
	r = vars[1];
	chi = vars[2];
	phi = vars[3];

	dt = vars[4];
	dr = vars[5];
	dchi = vars[6];
	dphi = vars[7];

	if (spherical == 0)
	{
		christoffel(r, chi, christ);
	}
	else
	{
		christoffel_spherical(r, chi, christ);
	}

	// 1st order diff eqs
	diffs[0] = dt;	 // dt
	diffs[1] = dr;	 // dr
	diffs[2] = dchi; // dchi
	diffs[3] = dphi; // dphi

	// 2nd order diff eqs for r and theta
	diffs[4] = -christ[0][0][0] * dt * dt - 2.0 * christ[0][0][1] * dt * dr - 2.0 * christ[0][0][2] * dt * dchi - 2.0 * christ[0][0][3] * dt * dphi - christ[0][1][1] * dr * dr - 2.0 * christ[0][1][2] * dr * dchi - 2.0 * christ[0][1][3] * dr * dphi - christ[0][2][2] * dchi * dchi - 2.0 * christ[0][2][3] * dchi * dphi - christ[0][3][3] * dphi * dphi; // d2t
	diffs[5] = -christ[1][0][0] * dt * dt - 2.0 * christ[1][0][1] * dt * dr - 2.0 * christ[1][0][2] * dt * dchi - 2.0 * christ[1][0][3] * dt * dphi - christ[1][1][1] * dr * dr - 2.0 * christ[1][1][2] * dr * dchi - 2.0 * christ[1][1][3] * dr * dphi - christ[1][2][2] * dchi * dchi - 2.0 * christ[1][2][3] * dchi * dphi - christ[1][3][3] * dphi * dphi; // d2r
	diffs[6] = -christ[2][0][0] * dt * dt - 2.0 * christ[2][0][1] * dt * dr - 2.0 * christ[2][0][2] * dt * dchi - 2.0 * christ[2][0][3] * dt * dphi - christ[2][1][1] * dr * dr - 2.0 * christ[2][1][2] * dr * dchi - 2.0 * christ[2][1][3] * dr * dphi - christ[2][2][2] * dchi * dchi - 2.0 * christ[2][2][3] * dchi * dphi - christ[2][3][3] * dphi * dphi; // d2chi
	diffs[7] = -christ[3][0][0] * dt * dt - 2.0 * christ[3][0][1] * dt * dr - 2.0 * christ[3][0][2] * dt * dchi - 2.0 * christ[3][0][3] * dt * dphi - christ[3][1][1] * dr * dr - 2.0 * christ[3][1][2] * dr * dchi - 2.0 * christ[3][1][3] * dr * dphi - christ[3][2][2] * dchi * dchi - 2.0 * christ[3][2][3] * dchi * dphi - christ[3][3][3] * dphi * dphi; // d2phi
}

void diffeqs_first(double vars[], double diffs[])
{
	double t, r, chi, phi;
	double dt, dr, dchi, dphi;
	double g[4][4], christ[4][4][4];
	double denom;

	t = vars[0];
	r = vars[1];
	chi = vars[2];
	phi = vars[3];

	dr = vars[4];
	dchi = vars[5];

	christoffel(r, chi, christ);
	metric(r, chi, g);

	denom = g[3][3]*g[0][0] - g[0][3]*g[0][3];
	dt = -(lz*g[0][3] + energy*g[3][3])/denom;
	dphi = (lz*g[0][0] + energy*g[0][3])/denom;

	/* 1st order diff eqs */
	diffs[0] = dt;	 // dt
	diffs[1] = dr;	 // dr
	diffs[2] = dchi; // dchi
	diffs[3] = dphi; // dphi

	/* 2nd order diff eqs for r and chi */
	diffs[4] = -christ[1][0][0] * dt * dt - 2.0 * christ[1][0][1] * dt * dr - 2.0 * christ[1][0][2] * dt * dchi - 2.0 * christ[1][0][3] * dt * dphi - christ[1][1][1] * dr * dr - 2.0 * christ[1][1][2] * dr * dchi - 2.0 * christ[1][1][3] * dr * dphi - christ[1][2][2] * dchi * dchi - 2.0 * christ[1][2][3] * dchi * dphi - christ[1][3][3] * dphi * dphi; // d2r
	diffs[5] = -christ[2][0][0] * dt * dt - 2.0 * christ[2][0][1] * dt * dr - 2.0 * christ[2][0][2] * dt * dchi - 2.0 * christ[2][0][3] * dt * dphi - christ[2][1][1] * dr * dr - 2.0 * christ[2][1][2] * dr * dchi - 2.0 * christ[2][1][3] * dr * dphi - christ[2][2][2] * dchi * dchi - 2.0 * christ[2][2][3] * dchi * dphi - christ[2][3][3] * dphi * dphi; // d2chi
}
