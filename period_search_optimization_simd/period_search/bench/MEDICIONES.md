# Mediciones de la ruta ASIMD

Registro de las optimizaciones evaluadas sobre `CalcStrategyAsimd`, con lo que
midio cada una y como aplicar las tres que sobreviven.

Fecha: 2026-08-21. Todas las cifras vienen de ejecutar los bancos de
[`bench/`](.) en el dispositivo real, no de estimaciones.

## Resumen

De las ocho propuestas de la lista original, **dos dan una mejora medible**. La
mayor de todas, sin embargo, no estaba en esa lista: es el bucle de derivadas
`Dg`, que resulto ser el 41% de `mrqcof`.

| Mejora | Efecto | Bit-exacta | Patch |
|---|---:|---|---|
| **9** — invertir el bucle `Dg` de `bright` | **−9,2%** de `bright` = −6,5% de `mrqcof` | si | [`patches/9-bucle-dg-invertido.patch`](patches/9-bucle-dg-invertido.patch) |
| **8** — izar punteros de fila en `mrqcof` | **−1,53%** de `mrqcof` | si | [`patches/8-izar-punteros-mrqcof.patch`](patches/8-izar-punteros-mrqcof.patch) |
| **1** — una sola division en `bright` | **−0,66%** de `bright` | si | [`patches/1-reciproco-bright.patch`](patches/1-reciproco-bright.patch) |
| ~~**1b**~~ — igual, sin refinamiento | ~~−4,58% de `bright`~~ | **descartada** | cambia el polo en 7 de 503 periodos, ver abajo |

**9 y 8 juntas dan −8,2% de `mrqcof`, sin cambiar un solo bit del resultado**
(`alpha`, `beta` y chi-cuadrado dan los mismos hashes). Las otras seis
propuestas no dieron nada, o salieron peor.

Medido sobre la **aplicacion completa**, 1+8+9 dan **−4,1%**. Ver
"Medida de extremo a extremo": el banco exagera la mejora 9, y el porque
importa.

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

| Variante | mediana | vs base | hash |
|---|---:|---:|---|
| baseline | 28,216 ms | — | `a3fb3240ff63555c` |
| 9: bucle `Dg` invertido | 26,369 ms | **−6,5%** | **igual** |
| **9 + 8** | **25,893 ms** | **−8,2%** | **igual** |
| 8: punteros de fila izados | 27,853 ms | −1,53% | igual |
| 5: `1/ave` izado | 28,241 ms | −0,6% (ruido) | distinto |

Los rangos de 9+8 (25,85-25,95 ms) no se acercan siquiera a los del baseline
(27,97-28,23 ms).

### Mejora 9 — el bucle `Dg` invertido

| Variante de `bright()` | mediana | vs base | hash `dyda` |
|---|---:|---:|---|
| baseline | 39,76 us | — | `7d254efad7444483` |
| **9: invertido, bit-exacto** | **36,12 us** | **−9,2%** | **igual** |
| invertido sin conservar el orden de reduccion | 35,43 us | −11,0% | distinto |

Renunciar a la exactitud aporta solo 1,8 puntos mas, asi que no compensa.

Verificada bit a bit en 14 configuraciones, barriendo `Lmax` de 2 a 8
(`Ncoef` = 9, 16, 25, 36, 49, 64, 81 — pares e impares, incluido el caso
limite `ncoef03 < 10` donde `cyklus1 = 0`) y `nrows` 4 y 6.

## Medida de extremo a extremo

Los numeros de arriba son de los bancos aislados. Ejecutando la **aplicacion
completa** sobre el `period_search_in` de referencia (6 periodos, `90 0.5 90.1`),
con el gobernador fijado en `performance` a 1,35 GHz y tres rondas alternas:

| Binario | mediana | rango | vs base |
|---|---:|---:|---:|
| base | 196 s | 195-197 | — |
| **1 + 8 + 9** (bit-exacto) | **188 s** | 187-192 | **−4,1%** |
| 1b + 8 + 9 | 182 s | 182-183 | −7,1% |

Los rangos no se solapan. **El fichero de salida es identico byte a byte en las
doce ejecuciones**, incluidas las de la variante 1b: mismo md5
`210a0526d6fe6c7964c774ee3da3864f`, mismos polos, misma chi-cuadrado impresa.

### El banco exagera la mejora 9

La ganancia real es bastante menor que la que predicen los bancos (−8,7% de
`mrqcof` para 1+8+9). Y las piezas no cuadran entre si:

- `exact` y `fast` se diferencian **solo** en el refinamiento de Markstein. Esos
  6 s son el 3,1% del total; en el banco esa misma diferencia es el 3,9% de
  `bright`. De ahi saldria que `bright` es el ~78% del tiempo de la aplicacion.
- Pero entonces la mejora 9 sola (−9,2% de `bright`) tendria que ahorrar ~14 s,
  y 1+8+9 juntas solo ahorran 8 s.

Es decir: **la mejora 1 rinde en la aplicacion lo que decia el banco, y la 9 no.**

La explicacion mas probable — sin verificar — es el estado de la cache.
`bench_bright` llama a `bright()` 20.480 veces seguidas con el mismo `Dg`, que
se queda caliente en L2. En la aplicacion, entre llamada y llamada se acumula
`alpha`/`beta` sobre una matriz de 54x54 y `curv()` reescribe `Dg` en cada
iteracion. Como **toda** la ganancia de la mejora 9 consiste en traer menos
lineas de cache, es justo la optimizacion mas sensible a esa diferencia.

Leccion general: un banco que mantiene el working set caliente **sobrestima
sistematicamente** una optimizacion de trafico de memoria. Para ese tipo de
cambio, la medida de extremo a extremo no es opcional.

### Proyeccion sobre un host real

Partiendo de las 17,7 h de media del host 814329:

| | tiempo | ahorro |
|---|---:|---:|
| Actual | 63.664 s (17,7 h) | — |
| **1 + 8 + 9** (bit-exacto) | 61.065 s (17,0 h) | **43 min** |
| 1b + 8 + 9 | 59.117 s (16,4 h) | 76 min |

Avisos: son 6 periodos de los ~550 de un workunit, medidos a frecuencia fija
mientras el host real corre con `ondemand`. El coste fijo de arranque entra
igual en las tres variantes, asi que si acaso **infravalora** la mejora relativa.

## El workunit completo descarta la 1b

Los 6 periodos de la prueba anterior daban salida identica en las tres
variantes. Ejecutando el **rango completo** del `period_search_in` de
referencia — 503 periodos, `base` y `1b+8+9` en paralelo, un nucleo cada uno —
la conclusion se invierte:

| | lineas |
|---|---:|
| identicas | 485 |
| solo el ultimo digito impreso (mismo polo) | 10 |
| chi-cuadrado distinta, mismo polo | 1 |
| **polo distinto** | **7** |

Siete de 503 intervalos (1,4%) reportan **un polo diferente**, con su periodo
ajustado y su chi-cuadrado. Uno de ellos no es un vecino:

```
linea 339:  base  P=96.47850756  chi2=5.600183  polo=( 85, 35)
            1b    P=96.47629620  chi2=5.614838  polo=(273, 53)
```

Esos 188 grados de diferencia en lambda son la ambiguedad del polo espejo. No
es un problema de precision: son minimos casi degenerados, y una perturbacion
de 1e-14 basta para decidir en cual cae el Levenberg-Marquardt. En unos casos
gana la 1b y en otros la base, sin patron.

**Conclusion: la 1b queda descartada.** Cambia el resultado cientifico, no solo
los bits. La 1 con refinamiento de Markstein da los mismos 0 ulp de siempre y
no tiene este problema.

Y la leccion de metodo: **6 periodos eran una muestra demasiado pequena**. Con
60 decisiones de polo la salida salia identica y la 1b parecia segura; hicieron
falta las 5.030 del rango completo para encontrar las siete que fallan.

### Bajo carga en los dos nucleos rinde mas

En esa ejecucion en paralelo, `1b+8+9` tardo 9.258 s y `base` 10.330 s: un
−10,4% frente al −7,1% medido en secuencial. Y esta **infravalorado**, porque
`base` termino 18 minutos despues y esos ultimos minutos los corrio solo.

Tiene sentido: con los dos nucleos compitiendo por la L2 compartida, una
optimizacion de trafico de memoria como la mejora 9 vale mas. Un host de BOINC
real corre con todos los nucleos ocupados, asi que la medida secuencial de un
solo nucleo probablemente **infravalora** lo que rinden estos parches en
produccion.

Medido despues, con el mismo protocolo — `base` y `1+8+9` en paralelo sobre los
503 periodos:

| Binario | tiempo | md5 de la salida |
|---|---:|---|
| base | 10.329 s | `d9a66743a6e9c7efb9c15243d30957df` |
| **1 + 8 + 9** | **9.594 s** (**−7,1%**) | `d9a66743a6e9c7efb9c15243d30957df` |

Confirma las dos cosas:

- **Bit-exacta a escala completa.** Los dos md5 coinciden entre si y ademas con
  el de la `base` de la ejecucion anterior, que fue independiente. Las 503
  lineas identicas, frente a las 7 que cambiaba la 1b.
- **Bajo carga de dos nucleos rinde −7,1%, no −4,1%.** La hipotesis de la
  contencion de L2 queda confirmada: la ganancia casi se duplica cuando los dos
  nucleos compiten por la cache compartida, que es como corre un host de BOINC
  real.

Y es una **cota inferior**: la variante rapida acabo 12 minutos antes, asi que
`base` corrio sola ese ultimo tramo. A igualdad de contencion la diferencia
seria algo mayor.

Confianza de la medida: `base` dio 10.330 s y 10.329 s en dos ejecuciones
independientes de casi tres horas. Un segundo de diferencia.

### Proyeccion revisada

Con −7,1%, que es la cifra aplicable a un host con todos los nucleos ocupados:

| | tiempo | ahorro | workunits/ano por nucleo |
|---|---:|---:|---:|
| Actual | 63.664 s (17,7 h) | — | 495 |
| **1 + 8 + 9** | 59.134 s (16,4 h) | **76 min** | **533** |

Es el mismo ahorro que daba la 1b, pero sin cambiar un solo bit del resultado.

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
| Bucle `Dg`: "no se arregla sin reestructurar `mrqcof`" | se arregla dentro de `bright`, −9,2% |

Tres hipotesis sobre por que el bucle `Dg` es lento (TLB, stride, latencia de
memoria) resultaron falsas las tres: ni compactar el stride ni el prefetch
cambiaron nada.

Ninguna estimacion basada en conteo de instrucciones y latencias publicadas
acerto. **En este nucleo, medir no es opcional.**

### Un falso positivo del 16%

Merece la pena dejarlo escrito. La primera version del bucle invertido midio
**−16,6%**, casi el doble de lo que rinde en realidad. Calculaba el numero de
vectores de `dyda` como `(ncoef0 - 3) >> 1`, y como `ncoef03 = 49` es **impar**,
eso trunca: cubria `dyda[0..47]` y **no calculaba `dyda[48]`**.

El hash lo detecto de inmediato, pero la diferencia se atribuyo al cambio de
orden de la reduccion, que era una explicacion plausible y falsa. La leccion es
que un hash distinto **no** se puede dar por explicado sin comprobarlo: solo un
hash *igual* prueba algo por si solo.

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

#### Cuanto se desvia la 1b, medido

Volcando los bits crudos de 512 geometrias (`PS_DUMP=1`) y comparando valor a
valor contra el baseline:

| | `ymod` | `dyda` |
|---|---|---|
| **1** (Markstein) | 0 ulp en 512/512 | **0 ulp en 29.184/29.184** |
| **1b** (sin refinar) | 1 ulp max., 2,7% de los valores | 45% difieren, hasta 2,5e-10 relativo |

Ese 2,5e-10 enganna. El peor caso cae sobre un componente de `dyda` que vale
4,5e-07, o sea **6e-11 veces el mayor del gradiente**: es ruido numerico sobre
una derivada que es cero a todos los efectos. Normalizado a la magnitud del
gradiente, que es lo que consume Levenberg-Marquardt:

```
error maximo sobre las 512 llamadas : 6,8e-14
mediana                             : 8,3e-16
```

Es decir: `ymod` correcto a 1 ulp y el gradiente a ~1e-14 relativo. La
chi-cuadrado converge al mismo punto fijo, porque ese punto lo define `ymod` y
no el camino de la iteracion.

**La decision es de quien la tome, no tecnica.** Para la ciencia no afecta.
Para la reproducibilidad si: es la primera divergencia respecto al binario
actual, y dos periodos con chi-cuadrado a menos de 1e-13 podrian intercambiarse
en el ranking. Las rutas AVX/SSE/CUDA ya difieren entre si, lo que **sugiere**
que el validador de Asteroids@home tolera esto por construccion — pero eso es
un argumento, no una comprobacion. Es lo que hay que mirar antes de desplegar
la 1b; la magnitud del error ya esta medida y es inofensiva.

Antes de aplicar cualquiera de las dos, `ps_div_exact` asume que no hay
subnormales en el residuo. Aqui `dnom = lmu + lmu0 > 2·TINY = 2e-8` y los
numeradores son de orden 1, asi que el margen es enorme.

### Comprobacion

```sh
cd period_search_optimization_simd/period_search/bench
make && scp -O bench_bright root@device:/tmp/ && ssh root@device '/tmp/bench_bright 40'
```

Con `PS_DUMP=1` el banco vuelca ademas los bits crudos de cada salida, lo que
permite comparar dos variantes valor a valor en vez de solo por el hash.

Con el patch sin modificar, `hash ymod` y `hash dyda` tienen que salir
**identicos** a los de antes de aplicarlo. Si cambian, la reestructuracion del
enmascarado esta mal. Con la variante rapida cambian por definicion.

## Como aplicar la mejora 9

Bit-exacta, −9,2% de `bright` = −6,5% de `mrqcof`.

```sh
cd <raiz del repo>
git apply period_search_optimization_simd/period_search/bench/patches/9-bucle-dg-invertido.patch
```

La operacion del bucle es

```
dyda[k] += suma_j  dbr[j] * Dg[fila_j][k]
```

un GEMV con intensidad aritmetica de 1 FMA por carga de 16 bytes. El bucle
original tiene `k` fuera y `j` dentro: para cada bloque de 10 coeficientes
barre las ~130 filas de `Dg` saltando `MAX_N_PAR+8` doubles entre una y otra, y
repite ese barrido seis veces. Cada fila se toca 11 veces a nivel de linea de
cache cuando solo tiene 7 lineas utiles.

Invertirlo — `j` fuera, `k` dentro — lee cada fila **de corrido y una sola
vez**, con los acumuladores de `dyda` en registros. Baja de 11 a 7 toques por
fila.

Que el coste lo fijan las lineas y no los FMA quedo demostrado por accidente:
tamanos de bloque que leen columnas de mas (10 y 16 vectores, que llegan hasta
el double 59 y 63 en lugar de parar en el 47) rinden **igual que el original**,
pese a ejecutar el mismo numero de FMA.

### Detalles que importan

- **`acc[]` tiene que quedarse en registros.** Eso exige que el tamano de
  bloque sea constante de compilacion, de ahi la descomposicion 20/8/4/2/1 en
  lugar de un bucle con tamano variable. Con 24 acumuladores GCC desborda a
  pila y el bucle pasa a ser un 49% **mas lento** que el original.
- **El prologo debe ser ligero.** Arrancar la reduccion sacando cuatro filas a
  la vez, o metiendo un condicional por acumulador, basta para que GCC deje de
  mantener `acc[]` en registros y la ganancia desaparece entera. La version
  final saca una sola multiplicacion.
- **La exactitud sale del orden de reduccion.** El bucle original suma sobre
  `j` en el orden 1,0,2,3,4,... por debajo de `cyklus1` y 0,1,2,3,4,... en la
  cola; es decir, solo intercambia las dos primeras. Se reproduce partiendo los
  bloques en esa frontera para que ninguno la cruce.

## Lo que queda sin atacar

- **El bucle de facetas y las llamadas a `matrix`/`phasec`**: el 29% de
  `mrqcof` que ningun banco ha medido todavia.
- **Procesar los puntos por lotes.** `bright()` se llama una vez por punto de
  datos **con el mismo `Dg`**, solo cambia `dbr`. Agrupar puntos convierte el
  GEMV en un GEMM y cada elemento de `Dg` serviria para muchos puntos, lo que
  atacaria lo que queda del bucle. Exige reestructurar `mrqcof` para calcular
  primero la geometria de todos los puntos y despues acumular, con riesgo real
  de cambiar resultados. Sin medir, y bastante mas caro que las mejoras de
  aqui.

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

`bench_bright [repeticiones] [nrows] [lmax]` — el tercer argumento permite
barrer `Ncoef = (lmax+1)^2`, que es como se comprobo que la mejora 9 es
bit-exacta para tamanos de problema distintos del de referencia.

Para comparar dos variantes, compilar las dos y alternarlas en la misma sesion
ssh. No comparar numeros tomados con horas de diferencia.
