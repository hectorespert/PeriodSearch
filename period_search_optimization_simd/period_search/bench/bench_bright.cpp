/* Banco aislado de CalcStrategyAsimd::bright().
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
 * Uso: bench_bright [repeticiones] [nrows] [lmax]
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
static int build_geometry(int nrows, int lmax)
{
    Lmax = Mmax = lmax;
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
    const long reps  = (argc > 1) ? atol(argv[1]) : 200;
    const int  nrows = (argc > 2) ? atoi(argv[2]) : 6;
    const int  lmax  = (argc > 3) ? atoi(argv[3]) : 6;

    build_geometry(nrows, lmax);

    /* ncoef03 == Ncoef  =>  ncoef = Ncoef + 5 + Nphpar */
    const int ncoef0 = Ncoef + 3;
    const int ncoef  = ncoef0 + 2 + Nphpar;

    std::vector<double> cg(ncoef + 2, 0.0);
    cg[1] = 0.10;                       /* tamano de la esfera */
    cg[ncoef0 - 2] = 0.55;              /* beta  (latitud del polo) */
    cg[ncoef0 - 1] = 1.10;              /* lambda */
    cg[ncoef0]     = 2.0 * PI * 24.0 / 90.0;   /* omega, periodo 90 h */
    cg[ncoef0 + 1] = 0.5;               /* parametros de fase */
    cg[ncoef0 + 2] = 0.1;
    cg[ncoef0 + 3] = -0.5;
    cg[ncoef - 1]  = log(0.5);          /* cl = exp(...) */
    cg[ncoef]      = 0.1;               /* cls */

    blmatrix(cg[ncoef0 - 2], cg[ncoef0 - 1]);

    auto strat = CreateAlignedShared<BenchStrategy>(64);

    /* Deformacion convexa suave. Los Legendre sin normalizar llegan a ~1e4, asi
       que se calibra la amplitud hasta que la razon de areas sea la de un
       asteroide real (factor < 3 entre la faceta mayor y la menor). */
    double amp = 1.0;
    for (int attempt = 0; attempt < 60; attempt++) {
        for (int i = 2; i <= Ncoef; i++) cg[i] = amp * sin(0.7 * i);
        strat->curv(cg, gl);
        double lo = 1e300, hi = 0.0;
        for (int i = 0; i < Numfac; i++) {
            const double ratio = gl.Area[i] / gl.Darea[i];   /* == exp(g) */
            if (ratio < lo) lo = ratio;
            if (ratio > hi) hi = ratio;
        }
        if (hi / lo < 3.0 && hi / lo > 1.2) { printf("amplitud=%.3g  Amax/Amin=%.2f\n", amp, hi / lo); break; }
        amp *= (hi / lo > 3.0) ? 0.5 : 1.6;
    }
    strat->curv(cg, gl);                /* rellena gl.Area y gl.Dg */

    /* Barrido de geometrias: rotacion + aspecto + angulo de fase variables */
    const int NS = 512;
    std::vector<double> ee(3 * NS), ee0(3 * NS), tt(NS);
    for (int s = 0; s < NS; s++) {
        const double u      = (double)s / (double)NS;
        const double aspect = 0.30 + 0.60 * u;
        const double phase  = 0.15 + 0.25 * u;
        ee0[3*s+0] = cos(aspect);            ee0[3*s+1] = 0.0;         ee0[3*s+2] = sin(aspect);
        ee [3*s+0] = cos(aspect)*cos(phase); ee [3*s+1] = sin(phase);  ee [3*s+2] = sin(aspect);
        double n = sqrt(ee[3*s]*ee[3*s] + ee[3*s+1]*ee[3*s+1] + ee[3*s+2]*ee[3*s+2]);
        ee[3*s] /= n; ee[3*s+1] /= n; ee[3*s+2] /= n;
        tt[s] = 2451000.0 + 0.37 * u;        /* barre varias rotaciones */
    }

    /* Fraccion de facetas visibles+iluminadas: comprueba que el reparto de
       ramas del bucle es realista (en un caso real ronda 0.2-0.4). */
    {
        long incl = 0, tot = 0;
        for (int s = 0; s < NS; s += 16) {
            for (int i = 0; i < Numfac; i++) {
                double lmu  = ee [3*s+0]*gl.Nor[0][i] + ee [3*s+1]*gl.Nor[1][i] + ee [3*s+2]*gl.Nor[2][i];
                double lmu0 = ee0[3*s+0]*gl.Nor[0][i] + ee0[3*s+1]*gl.Nor[1][i] + ee0[3*s+2]*gl.Nor[2][i];
                if (lmu > TINY && lmu0 > TINY) incl++;
                tot++;
            }
        }
        printf("Numfac=%d  Ncoef=%d  ncoef=%d  facetas incluidas=%.1f%%\n",
               Numfac, Ncoef, ncoef, 100.0 * (double)incl / (double)tot);
    }

    /* Calentamiento */
    for (int s = 0; s < NS; s++) {
        gl.xx1[1] = ee [3*s+0]; gl.xx1[2] = ee [3*s+1]; gl.xx1[3] = ee [3*s+2];
        gl.xx2[1] = ee0[3*s+0]; gl.xx2[2] = ee0[3*s+1]; gl.xx2[3] = ee0[3*s+2];
        strat->bright(tt[s], cg, ncoef, gl);
    }

    /* Medida: mediana de 5 pasadas */
    double best[5];
    uint64_t h_ymod = 0, h_dyda = 0;
    for (int run = 0; run < 5; run++) {
        h_ymod = 1469598103934665603ull; h_dyda = 1469598103934665603ull;
        const uint64_t t0 = ns_now();
        for (long r = 0; r < reps; r++) {
            for (int s = 0; s < NS; s++) {
                gl.xx1[1] = ee [3*s+0]; gl.xx1[2] = ee [3*s+1]; gl.xx1[3] = ee [3*s+2];
                gl.xx2[1] = ee0[3*s+0]; gl.xx2[2] = ee0[3*s+1]; gl.xx2[3] = ee0[3*s+2];
                strat->bright(tt[s], cg, ncoef, gl);
                hash_double(h_ymod, gl.ymod);
                for (int j = 0; j < ncoef; j++) hash_double(h_dyda, gl.dyda[j]);
            }
        }
        const uint64_t t1 = ns_now();
        best[run] = (double)(t1 - t0) / (double)(reps * NS);
    }
    std::sort(best, best + 5);

    printf("variante          : %s\n", kVariant);
    printf("llamadas          : %ld\n", reps * NS);
    printf("mediana           : %9.1f ns/llamada  (%.2f ns por par de facetas)\n",
           best[2], best[2] / (Numfac / 2.0));
    printf("min / max         : %9.1f / %.1f ns\n", best[0], best[4]);
    printf("hash ymod         : %016llx\n", (unsigned long long)h_ymod);
    printf("hash dyda         : %016llx\n", (unsigned long long)h_dyda);

    /* Con PS_DUMP=1 vuelca los bits crudos de cada salida, para comparar dos
       variantes valor a valor en lugar de solo por el hash. */
    if (getenv("PS_DUMP")) {
        for (int s = 0; s < NS; s++) {
            gl.xx1[1] = ee [3*s+0]; gl.xx1[2] = ee [3*s+1]; gl.xx1[3] = ee [3*s+2];
            gl.xx2[1] = ee0[3*s+0]; gl.xx2[2] = ee0[3*s+1]; gl.xx2[3] = ee0[3*s+2];
            strat->bright(tt[s], cg, ncoef, gl);
            uint64_t b; memcpy(&b, &gl.ymod, sizeof b);
            printf("y %016llx\n", (unsigned long long)b);
            for (int j = 0; j < ncoef; j++) {
                memcpy(&b, &gl.dyda[j], sizeof b);
                printf("d %016llx\n", (unsigned long long)b);
            }
        }
    }
    return 0;
}
