/* computes integrated brightness of all visible and illuminated areas
   and its derivatives

   8.11.2006 - Josef Durec
   25.3.2024 - Pavel Rosicky
*/

#include <math.h>
#include <cstdlib>
#include <cstdio>
#include <vector>
#include "globals.h"
#include "declarations.h"
#include "constants.h"
#include "CalcStrategyAsimd.hpp"
#include <arm_neon.h>

/* Cociente x/d correctamente redondeado a partir de r = fl(1/d) (Markstein).
   El residuo e = x - d*q0 es exacto porque el FMA no redondea el producto, y
   la correccion q0 + e*r da el mismo bit que FDIV. Cuesta 3 operaciones
   baratas en lugar de una division. */
static inline float64x2_t ps_div_exact(float64x2_t x, float64x2_t d, float64x2_t r)
{
    const float64x2_t q = vmulq_f64(x, r);
    const float64x2_t e = vfmsq_f64(x, d, q);
    return vfmaq_f64(q, e, r);
}

#define INNER_CALC \
	res_br = vaddq_f64(res_br, avx_pbr); \
	float64x2_t avx_sum1, avx_sum10, avx_sum2, avx_sum20, avx_sum3, avx_sum30; \
	\
    avx_sum1 = vmulq_f64(avx_Nor1, avx_de11); \
    avx_sum1 = vfmaq_f64(avx_sum1, avx_Nor2, avx_de21); \
    avx_sum1 = vfmaq_f64(avx_sum1, avx_Nor3, avx_de31); \
    \
    avx_sum10 = vmulq_f64(avx_Nor1, avx_de011); \
    avx_sum10 = vfmaq_f64(avx_sum10, avx_Nor2, avx_de021); \
    avx_sum10 = vfmaq_f64(avx_sum10, avx_Nor3, avx_de031); \
    \
    avx_sum2 = vmulq_f64(avx_Nor1, avx_de12); \
    avx_sum2 = vfmaq_f64(avx_sum2, avx_Nor2, avx_de22); \
    avx_sum2 = vfmaq_f64(avx_sum2, avx_Nor3, avx_de32); \
    \
    avx_sum20 = vmulq_f64(avx_Nor1, avx_de012); \
    avx_sum20 = vfmaq_f64(avx_sum20, avx_Nor2, avx_de022); \
    avx_sum20 = vfmaq_f64(avx_sum20, avx_Nor3, avx_de032); \
    \
    avx_sum3 = vmulq_f64(avx_Nor1, avx_de13); \
    avx_sum3 = vfmaq_f64(avx_sum3, avx_Nor2, avx_de23); \
    avx_sum3 = vfmaq_f64(avx_sum3, avx_Nor3, avx_de33); \
    \
    avx_sum30 = vmulq_f64(avx_Nor1, avx_de013); \
    avx_sum30 = vfmaq_f64(avx_sum30, avx_Nor2, avx_de023); \
    avx_sum30 = vfmaq_f64(avx_sum30, avx_Nor3, avx_de033); \
    \
    avx_sum1 = vmulq_f64(avx_sum1, avx_dsmu); \
    avx_sum2 = vmulq_f64(avx_sum2, avx_dsmu); \
    avx_sum3 = vmulq_f64(avx_sum3, avx_dsmu); \
    avx_sum10 = vmulq_f64(avx_sum10, avx_dsmu0); \
    avx_sum20 = vmulq_f64(avx_sum20, avx_dsmu0); \
    avx_sum30 = vmulq_f64(avx_sum30, avx_dsmu0); \
    \
    avx_dyda1 = vfmaq_f64(avx_dyda1, vaddq_f64(avx_sum1, avx_sum10), avx_Area); \
    avx_dyda2 = vfmaq_f64(avx_dyda2, vaddq_f64(avx_sum2, avx_sum20), avx_Area); \
    avx_dyda3 = vfmaq_f64(avx_dyda3, vaddq_f64(avx_sum3, avx_sum30), avx_Area); \
    \
    avx_d = vfmaq_f64(avx_d, avx_lmulmu0, avx_Area); \
    avx_d1 = vaddq_f64(avx_d1, avx_d1term);
// end of inner_calc

#define INNER_CALC_DSMU \
    avx_Area = vld1q_f64(&gl.Area[i]); \
    avx_dnom = vaddq_f64(avx_lmu, avx_lmu0); \
    avx_rdnom = vdivq_f64(avx_11, avx_dnom); \
    avx_lmulmu0 = vmulq_f64(avx_lmu, avx_lmu0); \
    avx_s = vmulq_f64(avx_lmulmu0, vaddq_f64(avx_cl, ps_div_exact(avx_cls, avx_dnom, avx_rdnom))); \
    avx_pdbr = vmulq_f64(vld1q_f64(&gl.Darea[i]), avx_s); \
    avx_pbr = vmulq_f64(avx_Area, avx_s); \
    avx_powdnom = ps_div_exact(avx_lmu0, avx_dnom, avx_rdnom); \
    avx_powdnom = vmulq_f64(avx_powdnom, avx_powdnom); \
    avx_dsmu = vfmaq_f64(vmulq_f64(avx_cls, avx_powdnom), avx_cl, avx_lmu0); \
    avx_powdnom = ps_div_exact(avx_lmu, avx_dnom, avx_rdnom); \
    avx_powdnom = vmulq_f64(avx_powdnom, avx_powdnom); \
    avx_dsmu0 = vfmaq_f64(vmulq_f64(avx_cls, avx_powdnom), avx_cl, avx_lmu); \
    avx_d1term = ps_div_exact(vmulq_f64(vmulq_f64(avx_Area, avx_lmu), avx_lmu0), avx_dnom, avx_rdnom);
// end of inner_calc_dsmu


#if defined __GNUG__ && !defined __clang__
__attribute__((__target__("arch=armv8-a+simd")))
#elif defined __GNUG__ && __clang__
// NOTE: The following generates warning: unsupported architecture 'armv8-a+simd' in the 'target' attribute string; 'target' attribute ignored [-Wignored-attributes]
// __attribute__((target("arch=armv8-a+simd")))
#endif

/**
 * @brief Computes integrated brightness of all visible and illuminated areas and its derivatives.
 *
 * This function calculates the integrated brightness of all visible and illuminated areas based on the provided time `t`,
 * coefficient vector `cg`, and global data. It also computes the derivatives of the brightness with respect to the coefficients.
 *
 * @param t The time at which the brightness is evaluated.
 * @param cg A reference to a vector of doubles containing the coefficients for the brightness calculation.
 * @param ncoef An integer representing the number of coefficients.
 * @param gl A reference to a globals structure containing necessary global data.
 *
 * @note The function modifies the global variables `ymod` and `dyda`.
 *
 * @date 8.11.2006
 * @author Josef Durec
 *
 * @date 25.3.2024 modified by Pavel Rosicky
 */
void CalcStrategyAsimd::bright(const double t, std::vector<double>& cg, const int ncoef, globals &gl)
{
   int i, j, k;
   incl_count = 0;
   double *ee = gl.xx1;
	double *ee0 = gl.xx2;

   ncoef0 = ncoef - 2 - Nphpar;
   cl = exp(cg[ncoef-1]); /* Lambert */
   cls = cg[ncoef];       /* Lommel-Seeliger */
   dot_product_new(ee, ee0, cos_alpha);
   alpha = acos(cos_alpha);
   for (i = 1; i <= Nphpar; i++)
      php[i] = cg[ncoef0+i];

   phasec(dphp,alpha,php); /* computes also Scale */

   matrix(cg[ncoef0],t,tmat,dtm);

   /* Directions (and derivatives) in the rotating system */

   for (i = 1; i <= 3; i++)
   {
      e[i] = 0;
      e0[i] = 0;
      for (j = 1; j <= 3; j++)
      {
         e[i] += tmat[i][j] * ee[j];
         e0[i] += tmat[i][j] * ee0[j];
         de[i][j] = 0;
         de0[i][j] = 0;
         for (k = 1; k <= 3; k++)
	 {
            de[i][j] += dtm[j][i][k] * ee[k];
            de0[i][j] += dtm[j][i][k] * ee0[k];
         }
      }
   }

   /*Integrated brightness (phase coefficients used later) */
   float64x2_t avx_e1 = vdupq_n_f64(e[1]);
   float64x2_t avx_e2 = vdupq_n_f64(e[2]);
   float64x2_t avx_e3 = vdupq_n_f64(e[3]);
   float64x2_t avx_e01 = vdupq_n_f64(e0[1]);
   float64x2_t avx_e02 = vdupq_n_f64(e0[2]);
   float64x2_t avx_e03 = vdupq_n_f64(e0[3]);
   float64x2_t avx_de11 = vdupq_n_f64(de[1][1]);
   float64x2_t avx_de12 = vdupq_n_f64(de[1][2]);
   float64x2_t avx_de13 = vdupq_n_f64(de[1][3]);
   float64x2_t avx_de21 = vdupq_n_f64(de[2][1]);
   float64x2_t avx_de22 = vdupq_n_f64(de[2][2]);
   float64x2_t avx_de23 = vdupq_n_f64(de[2][3]);
   float64x2_t avx_de31 = vdupq_n_f64(de[3][1]);
   float64x2_t avx_de32 = vdupq_n_f64(de[3][2]);
   float64x2_t avx_de33 = vdupq_n_f64(de[3][3]);
   float64x2_t avx_de011 = vdupq_n_f64(de0[1][1]);
   float64x2_t avx_de012 = vdupq_n_f64(de0[1][2]);
   float64x2_t avx_de013 = vdupq_n_f64(de0[1][3]);
   float64x2_t avx_de021 = vdupq_n_f64(de0[2][1]);
   float64x2_t avx_de022 = vdupq_n_f64(de0[2][2]);
   float64x2_t avx_de023 = vdupq_n_f64(de0[2][3]);
   float64x2_t avx_de031 = vdupq_n_f64(de0[3][1]);
   float64x2_t avx_de032 = vdupq_n_f64(de0[3][2]);
   float64x2_t avx_de033 = vdupq_n_f64(de0[3][3]);
   float64x2_t avx_Scale = vdupq_n_f64(Scale);

   float64x2_t avx_tiny = vdupq_n_f64(TINY);
   float64x2_t avx_cl = vdupq_n_f64(cl);
   float64x2_t avx_cl1 = vsetq_lane_f64(cl, vdupq_n_f64(1.0), 0);
   float64x2_t avx_cls = vdupq_n_f64(cls);
   float64x2_t avx_11 = vdupq_n_f64(1.0);
   float64x2_t res_br = vdupq_n_f64(0.0);
   float64x2_t avx_dyda1 = vdupq_n_f64(0.0);
   float64x2_t avx_dyda2 = vdupq_n_f64(0.0);
   float64x2_t avx_dyda3 = vdupq_n_f64(0.0);
   float64x2_t avx_d = vdupq_n_f64(0.0);
   float64x2_t avx_d1 = vdupq_n_f64(0.0);

   for (i = 0; i < Numfac; i += 2)
   {
      float64x2_t avx_lmu, avx_lmu0, cmpe, cmpe0, cmp;
      float64x2_t avx_Nor1 = vld1q_f64(&gl.Nor[0][i]);
      float64x2_t avx_Nor2 = vld1q_f64(&gl.Nor[1][i]);
      float64x2_t avx_Nor3 = vld1q_f64(&gl.Nor[2][i]);
      float64x2_t avx_s, avx_dnom, avx_dsmu, avx_dsmu0, avx_powdnom, avx_pdbr, avx_pbr;
      float64x2_t avx_rdnom, avx_lmulmu0, avx_d1term;
      float64x2_t avx_Area;

      avx_lmu = vmulq_f64(avx_e1, avx_Nor1);
      avx_lmu = vfmaq_f64(avx_lmu, avx_e2, avx_Nor2); \
      avx_lmu = vfmaq_f64(avx_lmu, avx_e3, avx_Nor3); \
      avx_lmu0 = vmulq_f64(avx_e01, avx_Nor1);
      avx_lmu0 = vfmaq_f64(avx_lmu0, avx_e02, avx_Nor2); \
      avx_lmu0 = vfmaq_f64(avx_lmu0, avx_e03, avx_Nor3); \

      cmpe = vreinterpretq_f64_u64(vcgtq_f64(avx_lmu, avx_tiny));
      cmpe0 = vreinterpretq_f64_u64(vcgtq_f64(avx_lmu0, avx_tiny));
      cmp = vreinterpretq_f64_u64(vandq_u64(vreinterpretq_u64_f64(cmpe), vreinterpretq_u64_f64(cmpe0)));
      int64x2_t cmp_int = vreinterpretq_s64_f64(cmp);
      int icmp = (vgetq_lane_s64(cmp_int, 0) & 1) | ((vgetq_lane_s64(cmp_int, 1) & 1) << 1);

	  if(icmp & 1)  //first and second or only first
      {
		 INNER_CALC_DSMU
		 if (icmp & 2) {
    		Dg_row[incl_count] = (float64x2_t*)&gl.Dg[i];

			float64_t tmp;
			vst1q_lane_f64(&tmp, avx_pdbr, 0);
			dbr[incl_count++] = vdupq_n_f64(tmp);

    		Dg_row[incl_count] = (float64x2_t*)&gl.Dg[i + 1];

         float64_t tmp2;
         vst1q_lane_f64(&tmp2, vextq_f64(avx_pdbr, avx_pdbr, 1), 0);
         dbr[incl_count++] = vdupq_n_f64(tmp2);
       } else {
         avx_pbr = vcombine_f64(vget_low_f64(avx_pbr), vdup_n_f64(0.0));
         avx_dsmu = vcombine_f64(vget_low_f64(avx_dsmu), vdup_n_f64(0.0));
         avx_dsmu0 = vcombine_f64(vget_low_f64(avx_dsmu0), vdup_n_f64(0.0));
         avx_lmulmu0 = vcombine_f64(vget_low_f64(avx_lmulmu0), vdup_n_f64(0.0));
         avx_d1term = vcombine_f64(vget_low_f64(avx_d1term), vdup_n_f64(0.0));

    		Dg_row[incl_count] = (float64x2_t*)&gl.Dg[i];

			float64_t tmp3;
			vst1q_lane_f64(&tmp3, avx_pdbr, 0);
         dbr[incl_count++] = vdupq_n_f64(tmp3);
		 }
		 INNER_CALC
	  }
	  else if (icmp & 2)
	  {
 		 INNER_CALC_DSMU

         avx_pbr = vcombine_f64(vget_high_f64(avx_pbr), vdup_n_f64(0.0));
         avx_dsmu = vcombine_f64(vdup_n_f64(0.0), vget_high_f64(avx_dsmu));
         avx_dsmu0 = vcombine_f64(vdup_n_f64(0.0), vget_high_f64(avx_dsmu0));
         avx_lmulmu0 = vcombine_f64(vdup_n_f64(0.0), vget_high_f64(avx_lmulmu0));
         avx_d1term = vcombine_f64(vdup_n_f64(0.0), vget_high_f64(avx_d1term));

         Dg_row[incl_count] = (float64x2_t*)&gl.Dg[i + 1];

         float64_t tmp4;
         vst1q_lane_f64(&tmp4, vextq_f64(avx_pdbr, avx_pdbr, 1), 0);
         dbr[incl_count++] = vdupq_n_f64(tmp4);

		 INNER_CALC
	  }
   }

   dbr[incl_count] = vdupq_n_f64(0.0);
   dbr[incl_count + 1] = vdupq_n_f64(0.0);
   dbr[incl_count + 2] = vdupq_n_f64(0.0);
   dbr[incl_count + 3] = vdupq_n_f64(0.0);
   Dg_row[incl_count] = Dg_row[0];
   Dg_row[incl_count + 1] = Dg_row[0];
   Dg_row[incl_count + 2] = Dg_row[0];
   Dg_row[incl_count + 3] = Dg_row[0];

   res_br = vpaddq_f64(res_br, res_br);
   vst1q_lane_f64(&gl.ymod, res_br, 0);

   /* Derivatives of brightness w.r.t. g-coefficients */
   /* Bucle invertido: j fuera (filas de Dg), k dentro (coeficientes).
      La version anterior recorria las incl_count filas de Dg una vez por cada
      bloque de coeficientes -- seis barridos, saltando MAX_N_PAR+8 doubles
      entre filas consecutivas. El coste de este bucle lo fijan las lineas de
      cache traidas de L2, no los FMA, asi que lo que importa es tocar cada
      fila las menos veces posible. Aqui cada fila se lee de corrido y los
      acumuladores de dyda viven en registros.

      El resultado es identico bit a bit al anterior: se reproduce su orden de
      reduccion sobre j, que suma 1,0,2,3,4,... para los coeficientes por
      debajo de cyklus1 y 0,1,2,3,4,... para los de la cola. Por eso los
      bloques se parten en esa frontera, para que ninguno la cruce. */
#define PS_DO_PRAGMA(x) _Pragma(#x)
#define PS_UNROLL_1(n)  PS_DO_PRAGMA(GCC unroll n)
#define PS_UNROLL(n)    PS_UNROLL_1(n)

/* Un bloque de KB vectores de dyda contra todas las filas de Dg. KB tiene que
   ser una constante de compilacion para que acc[] se quede en registros; de
   ahi que el recorrido lo descomponga en 20/8/4/2/1 en lugar de usar un
   tamano variable. */
#define PS_DG_BLOCK(KB)                                                       \
   {                                                                          \
      float64x2_t acc[KB];                                                    \
      const float64x2_t *rf = &Dg_row[jf][v0], *rs = &Dg_row[js][v0];         \
      const float64x2_t df = dbr[jf], ds = dbr[js];                           \
      int kk, jj;                                                             \
      PS_UNROLL(KB)                                                           \
      for (kk = 0; kk < KB; kk++) acc[kk] = vmulq_f64(df, rf[kk]);            \
      PS_UNROLL(KB)                                                           \
      for (kk = 0; kk < KB; kk++) acc[kk] = vfmaq_f64(acc[kk], ds, rs[kk]);   \
      for (jj = 2; jj < incl_count; jj++)                                     \
      {                                                                       \
         const float64x2_t pdbr = dbr[jj];                                    \
         const float64x2_t *row = &Dg_row[jj][v0];                            \
         PS_UNROLL(KB)                                                        \
         for (kk = 0; kk < KB; kk++)                                          \
            acc[kk] = vfmaq_f64(acc[kk], pdbr, row[kk]);                      \
      }                                                                       \
      PS_UNROLL(KB)                                                           \
      for (kk = 0; kk < KB; kk++)                                             \
         vst1q_f64(&gl.dyda[(v0 + kk) << 1], vmulq_f64(acc[kk], avx_Scale));  \
      v0 += KB;                                                               \
   }

#define PS_DG_RUN(END)                                                        \
   while (v0 < (END))                                                         \
   {                                                                          \
      const int left = (END) - v0;                                            \
      if      (left >= 20) PS_DG_BLOCK(20)                                    \
      else if (left >=  8) PS_DG_BLOCK(8)                                     \
      else if (left >=  4) PS_DG_BLOCK(4)                                     \
      else if (left >=  2) PS_DG_BLOCK(2)                                     \
      else                 PS_DG_BLOCK(1)                                     \
   }
   {
      const int ncoef03 = ncoef0 - 3;
      /* nv cubre los ncoef03 doubles de dyda; con ncoef03 impar hay que
         redondear hacia arriba. El double sobrante lo machaca despues el
         bloque de derivadas de rotacion, igual que hacia el bucle anterior. */
      const int nv  = (ncoef03 + 1) >> 1;
      const int nvh = ((ncoef03 / 10) * 10) >> 1;   /* frontera de orden */
      int v0 = 0;

      { const int jf = 1, js = 0; PS_DG_RUN(nvh) }   /* orden 1,0,2,3,... */
      { const int jf = 0, js = 1; PS_DG_RUN(nv)  }   /* orden 0,1,2,3,... */
   }

#undef PS_DG_RUN
#undef PS_DG_BLOCK
#undef PS_UNROLL
#undef PS_UNROLL_1
#undef PS_DO_PRAGMA

   /* Derivatives of brightness w.r.t. rotation parameters */
	avx_dyda1 = vpaddq_f64(avx_dyda1, avx_dyda2);
   avx_dyda1 = vmulq_f64(avx_dyda1, avx_Scale);
   vst1q_f64(&gl.dyda[ncoef0-3+1-1], avx_dyda1);  //unaligned memory because of odd index

   avx_dyda3 = vpaddq_f64(avx_dyda3, avx_dyda3);
   avx_dyda3 = vmulq_f64(avx_dyda3, avx_Scale);
   vst1q_f64(&gl.dyda[ncoef0-3+3-1], avx_dyda3); //unaligned memory because of odd index

   /* Derivatives of br. w.r.t. cl, cls */
   avx_d = vpaddq_f64(avx_d, avx_d1);
   avx_d = vmulq_f64(avx_d, avx_Scale);
   avx_d = vmulq_f64(avx_d, avx_cl1);
   vst1q_f64(&gl.dyda[ncoef-1-1], avx_d); //unaligned memory because of odd index

 /* Derivatives of br. w.r.t. phase function params. */
     for(i = 1; i <= Nphpar; i++)
       gl.dyda[ncoef0+i-1] = gl.ymod * dphp[i];

/*     dyda[ncoef0+1-1] = br * dphp[1];
     dyda[ncoef0+2-1] = br * dphp[2];
     dyda[ncoef0+3-1] = br * dphp[3];*/

   /* Scaled brightness */
	 gl.ymod *= Scale;
}
