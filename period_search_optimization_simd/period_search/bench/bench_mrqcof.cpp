/* Banco aislado de CalcStrategyAsimd::mrqcof().
 *
 * Reconstruye el estado global minimo que bright() necesita (geometria de la
 * esfera gaussiana, Nor/Darea/Area/Dg, coeficientes) usando las MISMAS rutinas
 * que la aplicacion (trifac, areanorm, sphfunc, curv), y luego llama a la
 * funcion real un numero configurable de veces.
 *
 * Sirve para dos cosas a la vez:
 *   - tiempo   : ns por llamada y por par de facetas
 *   - exactitud: checksum bit a bit de ymod y dyda[], para comparar variantes
 *
 * Uso: bench_mrqcof [repeticiones] [puntos] [nrows]
 */
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <cstdint>
#include <cstring>
#include <vector>
#include <algorithm>

#include "constants.h"
#include "arrayHelpers.hpp"
#include "globals.h"
#include "declarations.h"
#if defined(__aarch64__)
  #include "CalcStrategyAsimd.hpp"
  using BenchStrategy = CalcStrategyAsimd;
  static const char* kVariant = "ASIMD";
#else
  #include "CalcStrategyNone.hpp"
  using BenchStrategy = CalcStrategyNone;
  static const char* kVariant = "escalar (host)";
#endif
#include "SIMDHelpers.h"

/* --- globales que en produccion define period_search_BOINC.cpp --- */
int Lmax, Mmax, Niter, Lastcall, Ncoef, Numfac, Nphpar, Deallocate, n_iter;
double Ochisq, Chisq, Alamda, Alamda_incr, Alamda_start, Phi_0, Scale,
    Fc[MAX_N_FAC + 1][MAX_LM + 1], Fs[MAX_N_FAC + 1][MAX_LM + 1],
    Tc[MAX_N_FAC + 1][MAX_LM + 1], Ts[MAX_N_FAC + 1][MAX_LM + 1],
    Dsph[MAX_N_FAC + 1][MAX_N_PAR + 1],
    Blmat[4][4],
    Pleg[MAX_N_FAC + 1][MAX_LM + 1][MAX_LM + 1],
    Dblm[3][4][4];
std::vector<double> atry, beta, da;
CalcContext calcCtx;
SIMDSupport CPUopt;

globals gl;                        /* definicion del extern declarado en arrayHelpers.hpp */

/* FNV-1a sobre los bits crudos. Sumar los valores en un double no vale: el
   total crece hasta ~1e4 y su ulp acaba siendo mayor que las diferencias de
   1 ulp que se quieren detectar, asi que las absorbe. */
static inline void hash_double(uint64_t& h, double v)
{
    uint64_t bits;
    std::memcpy(&bits, &v, sizeof bits);
    h ^= bits;
    h *= 1099511628211ull;
}

static uint64_t ns_now()
{
    struct timespec t;
    clock_gettime(CLOCK_MONOTONIC, &t);
    return (uint64_t)t.tv_sec * 1000000000ull + t.tv_nsec;
}

/* Geometria identica a la de period_search_BOINC.cpp */
static int build_geometry(int nrows)
{
    Lmax = Mmax = 6;
    Nphpar = 3;
    Lastcall = 0;
    Phi_0 = 0.0;

    Ncoef = 0;
    for (int m = 0; m <= Mmax; m++)
        for (int l = m; l <= Lmax; l++) { Ncoef++; if (m != 0) Ncoef++; }

    std::vector<double> t(MAX_N_FAC + 1, 0.0), f(MAX_N_FAC + 1, 0.0);
    const double dth = PI / (2 * nrows);
    int k = 1;
    t[1] = 0; f[1] = 0;
    for (int i = 1; i <= nrows; i++) {
        double dph = PI / (2 * i);
        for (int j = 0; j <= 4 * i - 1; j++) { k++; t[k] = i * dth; f[k] = j * dph; }
    }
    for (int i = nrows - 1; i >= 1; i--) {
        double dph = PI / (2 * i);
        for (int j = 0; j <= 4 * i - 1; j++) { k++; t[k] = PI - i * dth; f[k] = j * dph; }
    }
    const int ndir = k + 1;
    t[ndir] = PI; f[ndir] = 0;

    Numfac = 8 * nrows * nrows;

    std::vector<std::vector<int>> ifp(MAX_N_FAC + 1, std::vector<int>(4, 0));
    trifac(nrows, ifp);

    std::vector<double> at(MAX_N_FAC + 1, 0.0), af(MAX_N_FAC + 1, 0.0);
    areanorm(t, f, ndir, Numfac, ifp, at, af, gl);
    sphfunc(Numfac, at, af);

    return ndir;
}


int main(int argc, char** argv)
{
    const long reps    = (argc > 1) ? atol(argv[1]) : 20;
    const int  npoints = (argc > 2) ? atoi(argv[2]) : 500;
    const int  nrows   = (argc > 3) ? atoi(argv[3]) : 6;

    build_geometry(nrows);

    const int ncoef0 = Ncoef + 3;
    const int ma     = Ncoef + 5 + Nphpar;      /* como en period_search_BOINC.cpp */

    /* Dos "curvas": una de luz relativa y la de regularizacion de convexidad */
    gl.Lcurves = 2;
    gl.Lpoints.assign(gl.Lcurves + 2, 0);
    gl.Inrel.assign(gl.Lcurves + 2, 0);
    gl.Lpoints[1] = npoints; gl.Inrel[1] = 1;   /* relativa: activa el bucle de 1/ave */
    gl.Lpoints[2] = 3;       gl.Inrel[2] = 0;   /* convexidad */

    const int ndata = npoints + 3;
    gl.maxLcPoints  = npoints;
    gl.maxDataPoints = ndata;
    gl.ytemp.assign(gl.maxLcPoints + 2, 0.0);
    gl.dytemp_sizeY = MAX_N_PAR + 1 + 4;
    gl.dytemp_sizeX = gl.maxLcPoints + 2;
    init_matrix(gl.dytemp, gl.dytemp_sizeX + 1, gl.dytemp_sizeY + 1, 0.0);
    gl.Weight.assign(ndata + 5, 1.0);
    gl.initializeVectors(MAX_N_PAR + 1, MAX_N_PAR + 8 + 1);

    std::vector<double> cg(ma + 2, 0.0);
    cg[1] = 0.10;
    cg[ncoef0 - 2] = 0.55;
    cg[ncoef0 - 1] = 1.10;
    cg[ncoef0]     = 2.0 * PI * 24.0 / 90.0;
    cg[ncoef0 + 1] = 0.5;
    cg[ncoef0 + 2] = 0.1;
    cg[ncoef0 + 3] = -0.5;
    cg[ma - 1]     = log(0.5);
    cg[ma]         = 0.1;

    blmatrix(cg[ncoef0 - 2], cg[ncoef0 - 1]);
    auto strat = CreateAlignedShared<BenchStrategy>(64);

    double amp = 1.0;
    for (int attempt = 0; attempt < 60; attempt++) {
        for (int i = 2; i <= Ncoef; i++) cg[i] = amp * sin(0.7 * i);
        strat->curv(cg, gl);
        double lo = 1e300, hi = 0.0;
        for (int i = 0; i < Numfac; i++) {
            const double r = gl.Area[i] / gl.Darea[i];
            if (r < lo) lo = r;
            if (r > hi) hi = r;
        }
        if (hi / lo < 3.0 && hi / lo > 1.2) { printf("amplitud=%.3g  Amax/Amin=%.2f\n", amp, hi / lo); break; }
        amp *= (hi / lo > 3.0) ? 0.5 : 1.6;
    }

    /* Datos: geometrias de observacion como en el banco de bright */
    std::vector<std::vector<double>> x1(ndata + 2, std::vector<double>(4, 0.0));
    std::vector<std::vector<double>> x2(ndata + 2, std::vector<double>(4, 0.0));
    std::vector<double> x3(ndata + 2, 0.0), y(ndata + 2, 1.0), sig(ndata + 2, 1.0);
    for (int p = 1; p <= ndata; p++) {
        const double u      = (double)(p - 1) / (double)ndata;
        const double aspect = 0.30 + 0.60 * u;
        const double phase  = 0.15 + 0.25 * u;
        double e0v[3] = { cos(aspect), 0.0, sin(aspect) };
        double ev [3] = { cos(aspect)*cos(phase), sin(phase), sin(aspect) };
        const double n = sqrt(ev[0]*ev[0] + ev[1]*ev[1] + ev[2]*ev[2]);
        for (int k = 0; k < 3; k++) { x1[p][k+1] = ev[k] / n; x2[p][k+1] = e0v[k]; }
        x3[p] = 2451000.0 + 0.37 * u;
        y[p]  = 1.0 + 0.2 * sin(11.0 * u);
    }
    for (int p = npoints + 1; p <= ndata; p++) sig[p] = 1.0 / 0.1;   /* peso de convexidad */

    /* ia igual que period_search_BOINC.cpp para una curva relativa */
    std::vector<int> ia(ma + 2, 0);
    for (int i = 0; i < Ncoef; i++) ia[i] = 1;
    ia[0] = 0;                                   /* fijo en ajuste relativo */
    ia[Ncoef + 0] = 1; ia[Ncoef + 1] = 1; ia[Ncoef + 2] = 1;      /* beta, lambda, periodo */
    ia[Ncoef + 3] = 1; ia[Ncoef + 4] = 1; ia[Ncoef + 5] = 1;      /* parametros de fase */
    ia[Ncoef + 3 + Nphpar + 1 - 1] = 0;          /* cl */
    ia[Ncoef + 3 + Nphpar + 2 - 1] = 0;          /* Lommel-Seeliger fijo */

    int mfit = 0, lastma = 0, lastone = 0;
    for (int j = 0; j < ma; j++) if (ia[j]) { mfit++; lastma = j; }
    for (int j = 1; j <= lastma; j++) { if (!ia[j]) break; lastone = j; }

    printf("Numfac=%d Ncoef=%d ma=%d mfit=%d lastone=%d lastma=%d  puntos=%d (%d + 3)\n",
           Numfac, Ncoef, ma, mfit, lastone, lastma, ndata, npoints);

    std::vector<double> bta(ma + 2, 0.0);
    double trial_chisq = 0.0;

    strat->mrqcof(x1, x2, x3, y, sig, cg, ia, ma, bta, mfit, lastone, lastma, trial_chisq, gl, false);
    printf("chi2 de referencia: %.17g\n", trial_chisq);

    double best[5];
    uint64_t h_beta = 0, h_alpha = 0;
    for (int run = 0; run < 5; run++) {
        h_beta = 1469598103934665603ull; h_alpha = 1469598103934665603ull;
        const uint64_t t0 = ns_now();
        for (long r = 0; r < reps; r++)
            strat->mrqcof(x1, x2, x3, y, sig, cg, ia, ma, bta, mfit, lastone, lastma, trial_chisq, gl, false);
        const uint64_t t1 = ns_now();
        best[run] = (double)(t1 - t0) / (double)reps;
        hash_double(h_beta, trial_chisq);
        for (int j = 0; j < mfit; j++) hash_double(h_beta, bta[j]);
        for (int j = 0; j < mfit; j++)
            for (int k = 0; k <= j; k++) hash_double(h_alpha, gl.alpha[j][k]);
    }
    std::sort(best, best + 5);

    printf("variante          : %s\n", kVariant);
    printf("llamadas          : %ld\n", reps);
    printf("mediana           : %12.1f ns/mrqcof  (%.1f ns por punto)\n", best[2], best[2] / ndata);
    printf("min / max         : %12.1f / %.1f ns\n", best[0], best[4]);
    printf("hash beta+chi2    : %016llx\n", (unsigned long long)h_beta);
    printf("hash alpha        : %016llx\n", (unsigned long long)h_alpha);
    return 0;
}
