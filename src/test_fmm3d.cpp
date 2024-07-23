#include <bla.hpp>
using namespace ngbla;

extern "C" {
void lfmm3d_s_c_g_(double *eps, int *nsource, double *source, double *charge,
                   double *pot, double *grad, int *ier);
}

double Kernel(Vec<3> px, Vec<3> py) {
  if (L2Norm2(px - py) == 0) return 0.0;
  return 1.0 / L2Norm(px - py);
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

  cout << "y = " << y << endl;
  cout << "grad = " << grad << endl;
  Apply (ptx, ptx, x, y);
  cout << "y ref = " << y << endl;
}
