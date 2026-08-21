# Banco de pruebas de la ruta ASIMD

Mide `CalcStrategyAsimd::bright()` y `mrqcof()` de forma aislada, sin ejecutar
un workunit completo.

> **Los resultados estan en [MEDICIONES.md](MEDICIONES.md)**: que optimizacion
> rindio cuanto, y como aplicar las tres que sobrevivieron. Los parches
> verificados estan en [`patches/`](patches/). La mayor es la 9, invertir el
> bucle de derivadas `Dg`: −8,2% de `mrqcof` junto con la 8, y bit-exacta. Un workunit tarda horas y su tiempo esta dominado por ruido; aqui un
ciclo de medida son segundos.

## Que hay

| Fichero | Que hace |
|---|---|
| `fdiv_probe.c` | Mide `vdivq_f64` frente a `vfmaq_f64` en el nucleo real. La **ratio** entre ambos es independiente de la frecuencia de reloj y es el factor que decide cuanto rinde eliminar divisiones del bucle de facetas. |
| `bench_bright.cpp` | Llama a la funcion `bright()` real (no una copia) sobre un estado reconstruido con las mismas rutinas que la aplicacion: `trifac`, `areanorm`, `sphfunc`, `curv`, `blmatrix`. |
| `bench_mrqcof.cpp` | Lo mismo un nivel mas arriba: monta una curva de luz relativa mas la de regularizacion de convexidad y llama a `mrqcof()`. Cubre lo que `bench_bright` no ve — la acumulacion de `alpha`/`beta`, el manejo de `dytemp` y el bucle de curvas relativas. |

`bench_bright` sirve para dos cosas a la vez:

- **Tiempo**: ns por llamada y por par de facetas, mediana de 5 pasadas.
- **Exactitud**: imprime un hash FNV-1a sobre los bits crudos de `ymod` y de
  `dyda[]`. Dos variantes que deban ser equivalentes bit a bit tienen que dar
  el mismo hash. Es una comprobacion mucho mas rapida y estricta que comparar
  la chi-cuadrado de un workunit, y recorre las cuatro ramas del bucle.

  El hash tiene que ser sobre los bits, no una suma en coma flotante. Acumular
  decenas de miles de valores de magnitud ~1 da un total de ~1e4, cuyo ulp ya
  es mayor que las diferencias de 1 ulp que se quieren detectar: una suma las
  absorbe y da falsos "identicos". Esto no es hipotetico, paso al comparar la
  variante sin ramas del bucle de facetas.

## Uso

```sh
# 1. libps_common.a tiene que existir
cd ../../../period_search_common && make -f Makefile_OpenWRT_aarch64

# 2. construir
cd ../period_search_optimization_simd/period_search/bench
make                                   # aarch64, contra ../../../../openwrt
make OPENWRT_DIR=/ruta/al/sdk          # o con otro SDK

# 3. ejecutar en el dispositivo (BusyBox no tiene sftp-server)
scp -O bench_bright fdiv_probe root@device:/tmp/
ssh root@device '/tmp/fdiv_probe && /tmp/bench_bright 200'
```

`bench_bright [repeticiones] [nrows] [lmax]` — por defecto 200, 6 y 6. Con
`nrows=6` y `lmax=6` sale `Numfac=288` y `Ncoef=49`, que es exactamente lo que
produce el `period_search_in` de referencia. `lmax` da `Ncoef = (lmax+1)^2` y
sirve para comprobar que una variante sigue siendo bit-exacta con tamanos de
problema distintos del de referencia, no solo con el.

`bench_mrqcof [repeticiones] [puntos] [nrows]` — por defecto 20, 500 y 6, que
reproduce `ma=57`, `mfit=54`, `lastone=lastma=54`. Ojo a esa igualdad: con una
curva relativa todos los parametros ajustables forman una tirada contigua de
unos, asi que el segundo bucle de `mrqcof` (el de "resto de parametros") no
llega a ejecutarse nunca.

## Validacion del propio arnes

`make host` construye la misma cosa para x86 usando la ruta escalar. No sirve
para medir SIMD, pero confirma que la geometria que levanta el arnes es la
correcta: `Numfac`, `Ncoef` y `ncoef` salen derivados, no puestos a mano, y
deben coincidir con los del workunit real.

La fraccion de facetas visibles **e** iluminadas que imprime (~46%) es el dato
que hay que vigilar: determina el reparto de ramas del bucle y cuantas facetas
entran en `incl_count`. Si esa fraccion se aleja mucho de la de un caso real, la
medida no representa nada.

La amplitud de la deformacion se autocalibra hasta que la razon entre la faceta
mayor y la menor queda por debajo de 3. Sin eso el arnes genera formas absurdas:
los polinomios de Legendre sin normalizar llegan a ~1e4 para l=6, asi que unos
coeficientes aparentemente pequenos disparan `exp(g)` a 1e80.

## Notas de medida

- El proceso BOINC de produccion compite por los nucleos. El `min / max` que
  imprime el banco hace visible ese ruido.
- El efecto que se busca al eliminar divisiones es de decenas de por ciento,
  muy por encima del ~1% de varianza entre pasadas.
