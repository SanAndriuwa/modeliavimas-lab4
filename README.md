# Modeliavimas-lab4

## Užduotis: skaitinis diferencijavimas ir integravimas

### 1. Analitiškai duotų funkcijų diferencijavimas

| Nr. | f(x) | x0 | f'(x0) | f''(x0) |
|---|---|---|---|---|
| 1 | x³-2x | 1 | 1 | 6 |
| 2 | sin(x) | π/3 | 0.5 | -√3/2 |
| 3 | exp(x) | 0 | 1 | 1 |

Naudokite `difapx()` ir apskaičiuokite absoliučias paklaidas, kai **h=0.1 ir h=0.01**:

$$f'(x_0)\approx\frac{f_1-f_{-1}}{2h},\qquad f'(x_0)\approx\frac{-f_2+8f_1-8f_{-1}+f_{-2}}{12h}.$$

$$f''(x_0)\approx\frac{f_1-2f_0+f_{-1}}{h^2},\qquad f''(x_0)\approx\frac{-f_2+16f_1-30f_0+16f_{-1}-f_{-2}}{12h^2}.$$

Čia f_j = f(x0+jh). PDF paskutinės formulės vardiklyje praleistas kvadratas; naudojama teisinga reikšmė **12h²**.

### 2. Lentelės duomenų diferencijavimas

| x1 | f1 | x2 | f2 | x3 | f3 |
|---|---|---|---|---|---|
| 0.8 | -1.0880 | 0.8472 | 0.7494 | -0.2 | 1.2214 |
| 0.9 | -1.0710 | 0.9472 | 0.8118 | -0.1 | 1.1052 |
| 1.0 | -1.0000 | 1.0472 | 0.8660 | 0.0 | 1.0000 |
| 1.1 | -0.8690 | 1.1472 | 0.9116 | 0.1 | 0.9048 |
| 1.2 | -0.6720 | 1.2472 | 0.9481 | 0.2 | 0.8187 |

Raskite pirmąją ir antrąją išvestines apie vidurinį tašką visomis keturiomis formulėmis. Interpoliuokite centrinius 3 ir visus 5 taškus Lagranžo daugianariais l2 ir l4, diferencijuokite juos ir palyginkite rezultatus.

### 3. Baigtinių integralų skaičiavimas

| Nr. | Integrandas | Ribos | Tikslioji reikšmė |
|---|---|---|---|
| 1 | x³-2x | [0, 2] | 0 |
| 2 | sin(x) | [0, π/2] | 1 |
| 3 | exp(-x) | [0, 1] | 1-exp(-1) |

Kai **N=4 ir N=8**, apskaičiuokite trapecijų, Simpsono ir Rombergo metodų paklaidas.

### 4. Netiesioginiai integralai

1. Integralas sin(x)/x nuo 0 iki ∞; baigtiniam skaičiavimui viršutinė riba 100.
2. Integralas exp(-x²) nuo 0 iki ∞; baigtiniam skaičiavimui viršutinė riba 100.

Raskite simbolines reikšmes ir palyginkite Simpsono, adaptyvųjį Simpsono, `quad()`, `quadl()` metodus, taip pat Gauss-Legendre pirmajam integralui ir Gauss-Hermite antrajam. Parametrai: **N=200, tol=10⁻⁴, IGL=20, a=0.001**.

Šaltinis: `Diff_Int.pdf`, `Int.zip` ir `difapx()` iš pateikto `pade.m`.

## Paleidimas

Atidarykite [lab4_main.m](lab4_main.m) ir paspauskite **Run**. Visi `.m` failai turi būti viename aplanke. `int()` daliai reikalingas **Symbolic Math Toolbox**. Be jo skriptas naudoja žinomas analitines reikšmes ir aiškiai praneša, kad simbolinis skaičiavimas nebuvo vykdytas.

Vykdomi visi variantai. Lentelėse atskirai pateikiamos kvadratūros paklaidos intervale [0.001,100] ir bendros paklaidos pradinio integralo [0,∞] atžvilgiu.

Paaiškinimas: [SOLUTION.md](SOLUTION.md).
