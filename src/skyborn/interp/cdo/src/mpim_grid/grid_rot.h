#ifndef GRID_ROT_H
#define GRID_ROT_H

double lamrot_to_lam(double phis, double rlas, double polphi, double pollam, double polgam);
double phirot_to_phi(double phis, double rlas, double polphi, double polgam);
void usvs_to_uv(double us, double vs, double phi, double rla, double polphi, double pollam, double *u, double *v);

#endif
