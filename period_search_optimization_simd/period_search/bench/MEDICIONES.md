# Mediciones de la ruta ASIMD

Registro de las optimizaciones evaluadas sobre `CalcStrategyAsimd`, con lo que
midio cada una y como aplicar las dos que sobreviven.

Fecha: 2026-08-21. Todas las cifras vienen de ejecutar los bancos de
[`bench/`](.) en el dispositivo real, no de estimaciones.

## Resumen

De ocho propuestas evaluadas, **dos dan una mejora medible**:

| Mejora | Efecto | Bit-exacta | Patch |
|---|---:|---|---|
| **8** — izar punteros de fila en `mrqcof` | **−1,53%** de `mrqcof` | si | [`patches/8-izar-punteros-mrqcof.patch`](patches/8-izar-punteros-mrqcof.patch) |
| **1** — una sola division en `bright` | **−0,66%** de `bright` | si | [`patches/1-reciproco-bright.patch`](patches/1-reciproco-bright.patch) |
| **1b** — igual, sin refinamiento | **−4,58%** de `bright` | no (~1 ulp) | variante de una linea sobre el patch 1 |

Las otras seis no dieron nada, o salieron peor. El detalle esta mas abajo.

Lo importante: **el objetivo real no estaba en la lista original**. El bucle de
derivadas `Dg` de `bright()` es el 41% de `mrqcof` y ninguna de las ocho lo
toca. Ver "Lo que queda".

## Plataforma y carga

- MediaTek MT7622, **2x Cortex-A53 a 1,1375 GHz** (in-order, 2 vias), OpenWrt 25.12.5, musl.
- GCC 14.3.0, `-march=armv8-a -O3 -std=c++2a`.
- Workunit de referencia (`period_search_in`): Lmax=Mmax=6 → `Ncoef=49`;
  nrows=6 → `Numfac=288`; 500 puntos de curva + 3 de convexidad;
  `ma=57`, `mfit=54`, `lastone=lastma=54`.

Costes medidos con `fdiv_probe`:

| Operacion | ns | ciclos aprox. |
|---|---:|---:|
| `vdivq_f64` (throughput) | 15,60 | ~18 |
| `vdivq_f64` (latencia) | 16,68 | ~19 |
| `vfmaq_f64` (throughput) | 2,23 | ~2,5 |
| **ratio div/fma** | | **7,0x** |

El divisor no esta segmentado (latencia ≈ throughput), pero **7x un FMA, no 25x**.
El manual del A53 sugiere cifras bastante peores; no coinciden con lo medido.

Reparto de ramas del bucle de facetas, sobre la geometria real:

```
pares: ninguno=49,4%  solo-carril0=5,2%  solo-carril1=4,6%  ambos=40,8%
```

Las facetas adyacentes comparten arista, asi que estan muy correlacionadas.
Solo el 9,8% de los pares son parciales.

## Metodologia

Tres herramientas, todas en este directorio:

1. **`fdiv_probe`** — coste de las instrucciones en el nucleo real. La ratio
   div/fma es independiente de la frecuencia.
2. **`bench_bright`** — llama a `bright()` real sobre estado reconstruido con
   las rutinas de la aplicacion.
3. **`bench_mrqcof`** — un nivel mas arriba; cubre el 30% de `mrqcof` que
   `bench_bright` no ve.

Dos tecnicas:

- **Ablacion**: para saber cuanto cuesta una parte, se desactiva y se mide la
  diferencia. Da la cota superior de lo que puede rendir optimizarla. Es la
  unica forma fiable de atribuir tiempo sin un profiler en el dispositivo.
- **Hash bit a bit**: FNV-1a sobre los bits crudos de las salidas. Determina si
  una variante es bit-exacta en segundos, sin ejecutar el workunit.

Cada comparacion se hace con **pasadas alternas** de las variantes (A, B, A, B,
...) para que la deriva termica o la carga del sistema afecte a las dos por
igual. Una diferencia solo se considera real si **los rangos no se solapan**.

El ruido tipico entre pasadas es de ±0,7%. Cualquier efecto por debajo de eso
es indistinguible.

## Donde esta el tiempo

Medido por ablacion:

```
mrqcof  28,285 ms  ....................................  100%
├─ bright  (500 llamadas x 39,55 us)  19,77 ms  ......   69,9%
│  ├─ bucle de derivadas Dg  ...........  23,0 ms/1000  58,5% de bright  = 40,9% de mrqcof
│  └─ resto (bucle de facetas, matrix, phasec)  ......   41,5% de bright = 29,0% de mrqcof
└─ resto  8,51 ms  ..................................    30,1%
   (acumulacion de alpha/beta, manejo de dytemp, bucle de curvas relativas)
```

Dos consecuencias:

- `bright` es el **69,9%** de `mrqcof`, no el ~90% que se venia suponiendo.
- El bucle `Dg` es lo mas caro con diferencia y no estaba en el radar.

Con `lastone == lastma` (caso de curva relativa), el segundo bucle de `mrqcof`
—el de "resto de parametros", lineas 195-227 y 266-297— **no se ejecuta nunca**.
Optimizar ahi no sirve de nada mientras todos los parametros ajustables formen
una tirada contigua de unos.

## Resultados

### `bright()` — 4 pasadas alternas, 20.480 llamadas cada una

| Variante | media | vs base | hash |
|---|---:|---:|---|
| baseline | 39.547 ns | — | `7d254efad7444483` |
| 1b: 1 division, sin refinamiento | 37.737 ns | **−4,58%** | distinto |
| *(ablacion: sin ninguna division)* | *37.312 ns* | *−5,5%* | *distinto* |
| 1: 1 division + Markstein | 39.285 ns | −0,66% | **igual** |
| 3: `vdupq_laneq_f64` | 39.647 ns | 0% | igual |
| 7: stride de `Dg` 208→72 doubles | 39.489 ns | 0% | igual |
| prefetch software en el bucle `Dg` | 39.952 ns | 0% / peor | igual |
| 2: `-mcpu=cortex-a53` | 39.348 ns | 0% | igual |
| 4: bucle de facetas sin ramas | 40.104 ns | **+1,0% peor** | distinto |

La variante 1b captura el 82% de lo que da la ablacion completa, que es el
techo teorico. Lo que queda es la division del reciproco, que no se puede quitar.

### `mrqcof()` — 3 pasadas alternas, 8 llamadas cada una

| Variante | media | vs base | hash |
|---|---:|---:|---|
| baseline | 28,285 ms | — | `a3fb3240ff63555c` |
| 8: punteros de fila izados | 27,853 ms | **−1,53%** | **igual** |
| 5: `1/ave` izado | 28,241 ms | −0,6% (ruido) | distinto |

La mejora 8 es la unica cuyos rangos **no se solapan** en ninguna pasada: la
peor ejecucion con punteros izados (27,908 ms) es mejor que la mejor sin ellos
(28,240 ms).

## Correcciones a las estimaciones previas

Merece la pena dejarlo escrito, porque el patron se repitio:

| Estimacion a priori | Medido |
|---|---|
| Divisiones: 60-84% del bucle de facetas | 5,5% de `bright` |
| `FDIV` vectorial: 30-60 ciclos | ~18 ciclos |
| Mejora 1: "1,6-2,5x" | −0,66% (bit-exacta) / −4,58% (1 ulp) |
| Mejora 2 (`-mcpu`): "5-10%" | 0% |
| Mejora 5: "4%, ~10% despues" | 0,6%, ruido |
| `bright` ~90% del total | 69,9% de `mrqcof` |
| Bucle `Dg`: prioridad secundaria | **58,5% de `bright`** |

Tres hipotesis sobre por que el bucle `Dg` es lento (TLB, stride, latencia de
memoria) resultaron falsas las tres: ni compactar el stride ni el prefetch
cambiaron nada.

Ninguna estimacion basada en conteo de instrucciones y latencias publicadas
acerto. **En este nucleo, medir no es opcional.**

## Como aplicar la mejora 8

Bit-exacta, −1,53% de `mrqcof`.

```sh
cd <raiz del repo>
git apply period_search_optimization_simd/period_search/bench/patches/8-izar-punteros-mrqcof.patch
```

Qué hace: en siete puntos de `mrqcof_asimd.cpp` saca el puntero de fila fuera
del bucle interno.

```c
// antes
for (m = 1; m <= l; m += 2) {
    float64x2_t avx_alpha = vld1q_f64(&alpha[j][k]);
    ...
    vst1q_f64(&alpha[j][k], avx_result);
}

// despues
double* arow = alpha[j].data();
for (m = 1; m <= l; m += 2) {
    float64x2_t avx_alpha = vld1q_f64(&arow[k]);
    ...
    vst1q_f64(&arow[k], avx_result);
}
```

Por que gana: `alpha` es un `vector<vector<double>>`, asi que cada acceso pasa
por dos indirecciones. Peor aun, el compilador no puede demostrar que
`alpha[j].data()` no aliasee con `gl.dyda`, de modo que **recarga el puntero de
fila despues de cada store**. Izarlo lo elimina. Lo mismo con `gl.dytemp[jp]`.

Los bucles de "resto de parametros" quedan sin tocar a proposito: con
`lastone == lastma` no se ejecutan. Si en el futuro alguna configuracion de `ia`
los activara, habria que izarlos tambien.

## Como aplicar la mejora 1

```sh
cd <raiz del repo>
git apply period_search_optimization_simd/period_search/bench/patches/1-reciproco-bright.patch
```

El patch deja `bright_asimd.cpp` con **una sola** `vdivq_f64` en lugar de
cuatro. Las tres divisiones de `INNER_CALC_DSMU` y la de `INNER_CALC` comparten
denominador, asi que se calcula `r = 1/dnom` una vez y se multiplica.

La cuarta parece irreducible porque usa el `dnom` ya enmascarado (las ramas
ponen `lmu=0`, `lmu0=1` en el carril muerto para evitar `0/0`). No lo es: se
calcula el termino con los valores **sin** enmascarar y se anula el *resultado*
en lugar de las entradas. El carril muerto aporta `+0.0` igual que antes, bit a
bit. Por eso el patch enmascara `avx_lmulmu0` y `avx_d1term` donde el original
enmascaraba `avx_lmu` y `avx_lmu0`.

### Las dos variantes

El patch tal cual usa el refinamiento de Markstein, que devuelve el cociente
**correctamente redondeado**, es decir identico bit a bit a `FDIV`:

```c
static inline float64x2_t ps_div_exact(float64x2_t x, float64x2_t d, float64x2_t r)
{
    const float64x2_t q = vmulq_f64(x, r);
    const float64x2_t e = vfmsq_f64(x, d, q);   /* residuo exacto gracias al FMA */
    return vfmaq_f64(q, e, r);
}
```

Cuesta 3 operaciones por cociente, y ahi se va casi toda la ganancia: **−0,66%**.

Para quedarse con el **−4,58%** hay que renunciar a la exactitud, sustituyendo
el cuerpo por:

```c
    (void)d;
    return vmulq_f64(x, r);   /* ~1 ulp de desviacion */
```

**La decision es de quien la tome, no tecnica.** El resultado se desvia en el
ultimo bit, lo que hace que la ruta ASIMD deje de coincidir con las demas. Las
rutas AVX/SSE/CUDA ya difieren entre si, asi que el validador de
Asteroids@home tolera esto por construccion — pero conviene comprobarlo antes
de desplegarlo, no despues.

Antes de aplicar cualquiera de las dos, `ps_div_exact` asume que no hay
subnormales en el residuo. Aqui `dnom = lmu + lmu0 > 2·TINY = 2e-8` y los
numeradores son de orden 1, asi que el margen es enorme.

### Comprobacion

```sh
cd period_search_optimization_simd/period_search/bench
make && scp -O bench_bright root@device:/tmp/ && ssh root@device '/tmp/bench_bright 40'
```

Con el patch sin modificar, `hash ymod` y `hash dyda` tienen que salir
**identicos** a los de antes de aplicarlo. Si cambian, la reestructuracion del
enmascarado esta mal. Con la variante rapida cambian por definicion.

## Lo que queda sin atacar

**El bucle de derivadas `Dg`** (`bright_asimd.cpp:264-352`), 41% de `mrqcof`.

No se arregla con micro-optimizacion, y hay evidencia: compactar el stride,
prefetch y `-mcpu` no movieron la aguja. El desensamblado son ~108
instrucciones por iteracion interna para 20-35 `fmla`, **sin spills** — GCC lo
compila bien.

El problema es estructural. La operacion es

```
dyda[k] += suma_j  dbr[j] * Dg[fila_j][k]
```

un GEMV donde **cada valor de `Dg` se carga y se usa exactamente una vez**:
intensidad aritmetica de 1 FMA por carga de 16 bytes. En un nucleo in-order de
2 vias eso esta limitado por emision de operaciones de memoria, y ninguna
reprogramacion del bucle lo cambia.

La reutilizacion existe, pero un nivel mas arriba: `bright()` se llama una vez
por punto de datos **con el mismo `Dg`**, solo cambia `dbr`. Procesar los puntos
por lotes convierte el GEMV en un GEMM y cada elemento de `Dg` sirve para muchos
puntos. Eso exige reestructurar `mrqcof` para calcular primero la geometria de
todos los puntos y despues acumular, con riesgo real de cambiar resultados.

Sin medir todavia: el 29% de `mrqcof` que es el bucle de facetas y las llamadas
a `matrix`/`phasec`.

## Reproducir

```sh
# 1. la libreria comun
cd period_search_common && make -f Makefile_OpenWRT_aarch64

# 2. los bancos
cd ../period_search_optimization_simd/period_search/bench
make                       # aarch64
make host                  # x86, ruta escalar, para validar el arnes

# 3. en el dispositivo
scp -O bench_bright bench_mrqcof fdiv_probe root@device:/tmp/
ssh root@device '/tmp/fdiv_probe; /tmp/bench_bright 40; /tmp/bench_mrqcof 8'
```

Para comparar dos variantes, compilar las dos y alternarlas en la misma sesion
ssh. No comparar numeros tomados con horas de diferencia.
