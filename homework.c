#include"fem.h"

# ifndef NOCONTACTITERATE

double femGrainsContactIterate(femGrains *myGrains, double dt, int iter)
{
	int i, j, iContact;
	int n = myGrains->n;


	double *x = myGrains->x;
	double *y = myGrains->y;
	double *m = myGrains->m;
	double *r = myGrains->r;
	double *vy = myGrains->vy;
	double *vx = myGrains->vx;
	double *dvBoundary = myGrains->dvBoundary;
	double *dvContacts = myGrains->dvContacts;
	double rIn = myGrains->radiusIn;
	double rOut = myGrains->radiusOut;

	double error = 0.0;
	double gap, rr, rx, ry, nx, ny, vn, dv, dvIn, dvOut;

	error = 0;
	iContact = 0;
	for (i = 0; i < n; i++) {
		for (j = i + 1; j < n; j++) {
			rx = (x[j] - x[i]);
			ry = (y[j] - y[i]);
			rr = sqrt(rx*rx + ry * ry);
			nx = rx / rr;
			ny = ry / rr;
			if (iter == 0) {
				dv = dvContacts[iContact];
			}
			else {
				vn = (vx[i] - vx[j])*nx + (vy[i] - vy[j])*ny;
				gap = rr - (r[i] + r[j]);
				dv = fmax(0.0, vn + dvContacts[iContact] - gap / dt);
				dv = dv - dvContacts[iContact];
				dvContacts[iContact] += dv;
				error = fmax(fabs(dv), error);
			}
			vx[i] -= dv * nx * m[j] / (m[i] + m[j]);
			vy[i] -= dv * ny * m[j] / (m[i] + m[j]);
			vx[j] += dv * nx * m[i] / (m[i] + m[j]);
			vy[j] += dv * ny * m[i] / (m[i] + m[j]);
			iContact++;
		}
	}

	for (i = 0; i < n; i++) {
		rr = sqrt(x[i] * x[i] + y[i] * y[i]);
		nx = x[i] / rr;
		ny = y[i] / rr;
		if (iter == 0) {
			dv = dvBoundary[i];
		}
		else {
			vn = vx[i] * nx + vy[i] * ny;
			gap = rOut - rr - r[i];
			dvOut = fmax(0.0, vn + dvBoundary[i] - gap / dt);
			gap = rr - rIn - r[i];
			dvIn = fmax(0.0, -vn - dvBoundary[i] - gap / dt);
			dv = dvOut - dvIn - dvBoundary[i];
			dvBoundary[i] += dv;
			error = fmax(fabs(dv), error);
		}
		vx[i] -= dv * nx;
		vy[i] -= dv * ny;
	}
	return error;
}

# endif
# ifndef NOUPDATE


void femGrainsUpdate(femGrains *myGrains, double dt, double tol, double iterMax)
{
	int i;
	int n = myGrains->n;

	double *x = myGrains->x;
	double *y = myGrains->y;
	double *m = myGrains->m;
	double *vy = myGrains->vy;
	double *vx = myGrains->vx;
	double gamma = myGrains->gamma;
	double gx = myGrains->gravity[0];
	double gy = myGrains->gravity[1];

	// 
	// -1- Calcul des nouvelles vitesses des grains sur base de la gravite et de la trainee
	//

	for (i = 0; i < n; i++) {
		double fx = m[i] * gx - gamma * vx[i];
		double fy = m[i] * gy - gamma * vy[i];
		vx[i] += fx * dt / m[i];
		vy[i] += fy * dt / m[i];
	}

	//
	// -2- Correction des vitesses pour tenir compte des contacts        
	//       

	int iter = 0;
	double error;

	do {
		error = femGrainsContactIterate(myGrains, dt, iter);
		iter++;
	} while ((error > tol / dt && iter < iterMax) || iter == 1);
	printf("iterations = %4d : error = %14.7e \n", iter - 1, error);

	//  
	// -3- Calcul des nouvelles positions sans penetrations de points entre eux
	//

	for (i = 0; i < n; ++i) {
		x[i] += vx[i] * dt;
		y[i] += vy[i] * dt;
	}
}

# endif