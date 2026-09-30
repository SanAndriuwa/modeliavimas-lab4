# Sprendimas: diferencijavimas ir integravimas

[lab4_main.m](lab4_main.m) apskaičiuoja visus užduoties variantus. `difapx.m` sprendžia nedidelę tiesinių lygčių sistemą finite-difference svoriams rasti; pirmosios išvestinės rezultatas dalijamas iš h, antrosios iš h².

## Diferencijavimas

Trijų taškų centrinės formulės paprastai turi O(h²) paklaidą, penkių taškų formulės O(h⁴). Kubinei funkcijai penkių taškų formulė yra tiksli iki apvalinimo paklaidos.

Pirmosios išvestinės absoliučios paklaidos, gautos MATLAB R2026a:

| Funkcija | h | 3 taškai | 5 taškai |
|---|---|---|---|
| x³-2x | 0.1 | 1.0000e-2 | 1.9984e-15 |
| x³-2x | 0.01 | 1.0000e-4 | 8.8818e-16 |
| sin(x) | 0.1 | 8.3292e-4 | 1.6647e-6 |
| sin(x) | 0.01 | 8.3333e-6 | 1.6667e-10 |
| exp(x) | 0.1 | 1.6675e-3 | 3.3373e-6 |
| exp(x) | 0.01 | 1.6667e-5 | 3.3337e-10 |

Lentelės duomenų išvestinės viduriniame taške:

| Duomenų rinkinys | f', 3 taškai | f', 5 taškai | f'', 3 taškai | f'', 5 taškai |
|---|---|---|---|---|
| 1 | 1.01000 | 1.00000 | 6.00000 | 6.00000 |
| 2 | 0.49900 | 0.49975 | -0.86000 | -0.8591667 |
| 3 | -1.00200 | -1.0004167 | 1.00000 | 0.9991667 |

Lagranžo l2 ir l4 daugianarių išvestinės sutampa su atitinkamomis baigtinių skirtumų reikšmėmis iki apvalinimo. l2 naudojami centriniai 3 taškai; l4 naudojami visi 5 taškai.

PDF 3 rinkinio lentelė atitinka **exp(-x)**, nors pirmojoje analitinių funkcijų užduotyje nurodyta **exp(x)**. Abi dalys skaičiuojamos pagal savo pateiktus duomenis. Suapvalinti lentelių duomenys riboja išvestinių tikslumą.

## Baigtiniai integralai

Trapecijų metodas sumuoja trapecijų plotus. Simpsono metodas naudoja svorius 1,4,2,...,4,1. Rombergo lentelė patikslina trapecijų rezultatą Richardsono ekstrapoliacija; paskutinė jos eilutė naudoja tiksliai N intervalų.

| Integralas | N | Trapecijų paklaida | Simpsono paklaida | Rombergo paklaida |
|---|---|---|---|---|
| x³-2x | 4 | 0.25 | 0 | 0 |
| x³-2x | 8 | 0.0625 | 0 | 0 |
| sin(x) | 4 | 1.28842e-2 | 1.34585e-4 | 8.43453e-6 |
| sin(x) | 8 | 3.21483e-3 | 8.29552e-6 | 8.14402e-9 |
| exp(-x) | 4 | 3.28887e-3 | 1.36165e-5 | 3.16180e-7 |
| exp(-x) | 8 | 8.22859e-4 | 8.55776e-7 | 1.23418e-10 |

## Integralai iki begalybės

Analitinės reikšmės:

$$\int_0^\infty\frac{\sin x}{x}\,dx=\frac\pi2,\qquad\int_0^\infty e^{-x^2}\,dx=\frac{\sqrt\pi}{2}.$$

Skaitiniai metodai skaičiuoja intervale **[0.001,100]**. Jo tikrosios reikšmės yra 1.561225466945 ir 0.885226925786. Skirtumas nuo integralo [0,∞] yra atitinkamai apie **0.00957086** ir **0.001**. Ši ribų pakeitimo paklaida atskirta nuo kvadratūros paklaidos.

Kvadratūros absoliučios paklaidos MATLAB R2026a, laikantis užduoties parametrų:

| Metodas | sin(x)/x | exp(-x²) |
|---|---|---|
| Simpsonas, N=200 | 3.20996e-6 | 3.69291e-5 |
| Adaptyvus Simpsonas | 3.76691e-8 | 6.30259e-8 |
| quad | 2.00064e-1 | 1.29146e-6 |
| quadl | 6.40002e-9 | 4.26506e-9 |
| Gauss-Legendre / Gauss-Hermite, 20 mazgų | 8.27685e-2 | 1.11022e-16 |

Adaptyvus Simpsonas lygina viso intervalo įvertį su dviejų pusių suma, dalija skirtumą iš 15 ir rekursiškai dalija intervalą. Pradžioje [a,b] suskaidomas į ne ilgesnes kaip 1 dalis: kitaip retai išdėstyti pirmieji taškai gali nepastebėti sin(x)/x svyravimų ir sustabdyti algoritmą per anksti. Bendras tol paskirstomas pradinėms dalims ir rekursijos šakoms.

`quad()` tiesioginis iškvietimas šiame svyruojančiame pavyzdyje davė didelę paklaidą, nors tol=10⁻⁴. Paklaidos įvertis nėra garantija, jei pradiniai taškai praleidžia svyravimus. Šis rezultatas paliktas palyginimo lentelėje. 20 Gauss-Legendre mazgų visame ilgame intervale taip pat nepakanka geram rezultatui.

Gauss-Hermite svoriai skirti **exp(-x²)** svoriui visoje realioje ašyje. Antrojo integralo pusės ašies reikšmė gaunama iš pusės svorių sumos; atimamas intervalas [0,a] ir uodega [b,∞]. Tai ištaiso pateikto `Gauss_Hermite.m` netinkamą baigtinio intervalo pakeitimą. Gauss metodų mazgai ir svoriai randami iš simetrinių rekurencijos matricų.

## Patikrinimas

Visi variantai paleisti MATLAB R2026a. Adaptyviojo metodo paklaida abiem integralams mažesnė už 10⁻⁴. Statinis analizatorius pateikia tik dvi rekomendacijas pakeisti senas `quad()` ir `quadl()` funkcijas į `integral()`; jos paliktos, nes aiškiai reikalaujamos užduotyje.

Symbolic Math Toolbox šiame kompiuteryje nėra, todėl **`int()` faktiškai nebuvo vykdyta**. Skripte jos iškvietimas išsaugotas ir automatiškai naudojamas, kai toolbox yra. Be jo skaičiavimui naudojamos pateiktos analitinės reikšmės; sin(x)/x baigtinio intervalo etalonas skaičiuotas `integral()` su 10⁻¹² tolerancijomis.
