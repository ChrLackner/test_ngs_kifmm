
#include <bla.hpp>
using namespace ngbla;

extern "C" {
void lfmm3d_s_c_g_(double *eps, int *nsource, double *source, double *charge,
                   double *pot, double *grad, int *ier);
}

double Kernel(Vec<3> px, Vec<3> py) {
  if (L2Norm2(px - py) == 0) return 0.0;
  return 1.0 / (4 * M_PI) / L2Norm(px - py);
}

void Apply(const Array<Vec<3>> &ptx, const Array<Vec<3>> &pty, FlatVector<> x,
           FlatVector<> y) {
  y = 0;
  for (size_t ix = 0; ix < ptx.Size(); ix++)
    for (size_t iy = 0; iy < pty.Size(); iy++)
      y(iy) += Kernel(ptx[ix], pty[iy]) * x(ix);
}

int main() {
  int n = 10;
  Array<Vec<3>> ptx(n);

  for (size_t i = 0; i < ptx.Size(); i++)
    for (size_t j = 0; j < 3; j++) ptx[i](j) = double(rand()) / RAND_MAX;

  Vector<> x(n), y(n);
  Vector<Vec<3>> grad(n);
  grad = 0;

  for (size_t i = 0; i < x.Size(); i++) x(i) = double(rand()) / RAND_MAX;
  y = 0;

  Vector<uintptr_t> expansion_order(1);
  expansion_order = 5;  // expansion order of FMM

  double eps = 0.5e-6;
  int ier = 0;

  lfmm3d_s_c_g_(&eps, &n, ptx[0].Data(), x.Data(), y.Data(), grad[0].Data(),
                &ier);

  // Apply (ptx, ptx, x, y);
  cout << "y = " << y << endl;
}

// #include "math.h"
// #include "stdio.h"
// #include "stdlib.h"

// extern "C" {
// #include "complex.h"
// #include "cprini.h"
// #include "lfmm3d_c.h"
// }

// int main(int argc, char **argv) {
//   // cprin_init("stdout", "fort.13");
//   // cprin_skipline(2);

//   int ns = 2000;
//   int ier = 0;
//   double *source = (double *)malloc(3 * ns * sizeof(double));

//   double *charge = (double *)malloc(ns * sizeof(double));

//   double *pot = (double *)malloc(ns * sizeof(double));
//   double *grad = (double *)malloc(3 * ns * sizeof(double));

//   // initialize arrays
//   for (int i = 0; i < ns; i++) {
//     source[3 * i] = pow(rand01(), 2);
//     source[3 * i + 1] = pow(rand01(), 2);
//     source[3 * i + 2] = pow(rand01(), 2);

//     charge[i] = rand01();
//   }

//   double eps = 0.5e-6;
//   // cprin_message("this code is an example c driver.");
//   // cprin_message("on output, the code prints sample pot,grad");
//   // cprin_skipline(2);

//   // call the fmm routine
//   lfmm3d_s_c_g_(&eps, &ns, source, charge, pot, grad, &ier);

//   // cprind("pot=", pot, 12);
//   // cprind("grad=", grad, 12);

//   return 0;
// }
