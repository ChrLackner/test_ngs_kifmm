
#include <bla.hpp>
using namespace ngbla;

#include <kifmm_rs.h>

double Kernel(Vec<3> px, Vec<3> py) {
  if (L2Norm2(px - py) == 0) return 0.0;
  return 1.0 / (4 * M_PI) / L2Norm(px - py);
}

void Apply(const Array<Vec<3>>& ptx, const Array<Vec<3>>& pty, FlatVector<> x,
           FlatVector<> y) {
  y = 0;
  for (size_t ix = 0; ix < ptx.Size(); ix++)
    for (size_t iy = 0; iy < pty.Size(); iy++)
      y(iy) += Kernel(ptx[ix], pty[iy]) * x(ix);
}

int main() {
  int n = 1e6;
  Array<Vec<3>> ptx(n);

  for (size_t i = 0; i < ptx.Size(); i++)
    for (size_t j = 0; j < 3; j++) ptx[i](j) = double(rand()) / RAND_MAX;

  Vector<> x(n), y(n);

  for (size_t i = 0; i < x.Size(); i++) x(i) = double(rand()) / RAND_MAX;
  y = 0;

  Vector<uintptr_t> expansion_order(1);
  expansion_order = 5;  // expansion order of FMM

  static Timer t1("KiFMM setup");
  t1.Start();
  // auto fmm = laplace_blas_f64(
  auto fmm = laplace_fft_f64(
      expansion_order.Data(), expansion_order.Size(), ptx.Data()->Data(),
      ptx.Size() * 3, ptx.Data()->Data(), ptx.Size() * 3, x.Data(), x.Size(),
      true,  // prune empty
      150,   // n_crit
      0,     // depth
      1);
  // 1e-4,  // singular_value_threshold
  // false, 0, 0, 0);
  t1.Stop();
  static Timer t2("KiFMM apply");
  t2.Start();
  // evaluate_laplace_blas_f64(fmm, true);
  evaluate_laplace_fft_f64(fmm, true);
  t2.Stop();

  cout << "kiFMM setup time: " << t1.GetTime() << endl;
  cout << "kiFMM apply time: " << t2.GetTime() << endl;
  cout << "kiFMM total time: " << t1.GetTime() + t2.GetTime() << endl;

  // Apply (ptx, ptx, x, y);
  // cout << "y = " << y << endl;
}
