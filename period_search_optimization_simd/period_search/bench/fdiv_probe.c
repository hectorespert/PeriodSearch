/* Mide el coste real de vdivq_f64 frente a vfmaq_f64 en el nucleo donde se ejecuta.
 * Lo que importa es la RATIO: es independiente de la frecuencia del reloj, y es
 * exactamente el factor que decide cuanto rinde la mejora 1 de bright_asimd.cpp.
 */
#include <arm_neon.h>
#include <stdint.h>
#include <stdio.h>
#include <time.h>

static uint64_t ns_now(void)
{
    struct timespec t;
    clock_gettime(CLOCK_MONOTONIC, &t);
    return (uint64_t)t.tv_sec * 1000000000ull + t.tv_nsec;
}

#define REPS 1000000L

static volatile double sink;

/* Throughput: 4 operaciones independientes por iteracion. */
static double thr_div(void)
{
    float64x2_t k = vdupq_n_f64(1.0), one = vdupq_n_f64(1.0);
    float64x2_t d0 = vdupq_n_f64(3.0), d1 = vdupq_n_f64(5.0);
    float64x2_t d2 = vdupq_n_f64(7.0), d3 = vdupq_n_f64(11.0);
    float64x2_t s0 = vdupq_n_f64(0), s1 = vdupq_n_f64(0);
    float64x2_t s2 = vdupq_n_f64(0), s3 = vdupq_n_f64(0);
    uint64_t t0 = ns_now();
    for (long i = 0; i < REPS; i++) {
        k  = vaddq_f64(k, one);
        s0 = vaddq_f64(s0, vdivq_f64(k, d0));
        s1 = vaddq_f64(s1, vdivq_f64(k, d1));
        s2 = vaddq_f64(s2, vdivq_f64(k, d2));
        s3 = vaddq_f64(s3, vdivq_f64(k, d3));
    }
    uint64_t t1 = ns_now();
    sink = vgetq_lane_f64(vaddq_f64(vaddq_f64(s0, s1), vaddq_f64(s2, s3)), 0);
    return (double)(t1 - t0) / (double)(REPS * 4);
}

static double thr_fma(void)
{
    float64x2_t k = vdupq_n_f64(1.0), one = vdupq_n_f64(1e-9);
    float64x2_t d0 = vdupq_n_f64(3.0), d1 = vdupq_n_f64(5.0);
    float64x2_t d2 = vdupq_n_f64(7.0), d3 = vdupq_n_f64(11.0);
    float64x2_t s0 = vdupq_n_f64(0), s1 = vdupq_n_f64(0);
    float64x2_t s2 = vdupq_n_f64(0), s3 = vdupq_n_f64(0);
    uint64_t t0 = ns_now();
    for (long i = 0; i < REPS; i++) {
        k  = vaddq_f64(k, one);
        s0 = vfmaq_f64(s0, k, d0);
        s1 = vfmaq_f64(s1, k, d1);
        s2 = vfmaq_f64(s2, k, d2);
        s3 = vfmaq_f64(s3, k, d3);
    }
    uint64_t t1 = ns_now();
    sink = vgetq_lane_f64(vaddq_f64(vaddq_f64(s0, s1), vaddq_f64(s2, s3)), 0);
    return (double)(t1 - t0) / (double)(REPS * 4);
}

/* Latencia: una sola cadena dependiente. */
static double lat_div(void)
{
    float64x2_t x = vdupq_n_f64(1.234567), c = vdupq_n_f64(1.000001);
    uint64_t t0 = ns_now();
    for (long i = 0; i < REPS; i++) x = vdivq_f64(c, x);
    uint64_t t1 = ns_now();
    sink = vgetq_lane_f64(x, 0);
    return (double)(t1 - t0) / (double)REPS;
}

int main(void)
{
    thr_div(); thr_fma(); lat_div();            /* calentamiento */

    double D = thr_div(), F = thr_fma(), L = lat_div();

    printf("vdivq_f64  throughput : %7.3f ns/op\n", D);
    printf("vdivq_f64  latencia   : %7.3f ns/op\n", L);
    printf("vfmaq_f64  throughput : %7.3f ns/op\n", F);
    printf("ratio div/fma         : %7.2f x\n\n", D / F);

    /* Cuerpo del bucle de facetas de bright(): ~60 ops baratas + N divisiones. */
    double viejo = 4.0 * D + 60.0 * F;   /* 4 vdivq */
    double nuevo = 1.0 * D + 72.0 * F;   /* 1 vdivq + refinamiento Markstein */
    printf("Estimacion para el bucle de bright():\n");
    printf("  actual (4 div)      : %7.1f ns/par-de-facetas\n", viejo);
    printf("  con 1 div + Markstein: %6.1f ns/par-de-facetas\n", nuevo);
    printf("  speedup del bucle   : %7.2f x\n", viejo / nuevo);
    return 0;
}
