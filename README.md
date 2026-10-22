# Modeliavimas-lab4

Atnaujinta pagal `4. Skaitinio diferencijavimo ir integravimo metodų taikymas 2022.pdf` (2026-10-22).

Atidarykite `lab4_main.m` ir paspauskite Run. Visi `.m` failai turi būti viename aplanke. Simboliniam integralui reikalingas Symbolic Math Toolbox.

Be šio toolbox programa veikia su žinoma π/2 reikšme ir aiškiai praneša, kad int nebuvo vykdytas.

Skripto pradžioje nustatykite `variantas = 1`, 2 arba 3 pagal dėstytojo nurodymą. Numatyta 1, nes PDF asmeninio varianto formulės nepateikia.

Programa pateikia užduoties lenteles:

1. Pirmosios ir antrosios išvestinių absoliučios paklaidos, naudojant 3 ir 5 taškus, h=0.1 ir 0.01.
2. Lentelės duomenų išvestinės ties viduriniu tašku bei palyginimas su 2 ir 4 laipsnio Lagranžo daugianarių išvestinėmis.
3. Baigtinio integralo trapecijų, Simpsono ir Rombergo paklaidos, kai N=4 ir 8.
4. Vienintelis naujos užduoties netiesioginis integralas sin(x)/x: simbolinis rezultatas ir Simpsono, adaptyviojo Simpsono, quad, quadl, Gauss-Legendre palyginimas.

Naudojami N=200, tol=1e-4, IGL=20, a=0.001 ir b=100. Pateikiama tiek bendra paklaida iki tikslaus begalinio integralo, tiek skaitinė paklaida baigtiniame intervale. Šiam darbui grafikai užduotyje neprašomi.

Paaiškinimas: [SOLUTION.md](SOLUTION.md).
