# Sprendimo paaiškinimas

## Diferencijavimas

Vietinė funkcija `isvestines()` tiesiogiai užrašo keturias PDF formules. Paskutinės formulės vardiklis yra 12h². Rezultatų tvarka: pirmoji 3 taškų, pirmoji 5 taškų, antroji 3 taškų, antroji 5 taškų išvestinė.

Analitinės funkcijos: x³−2x ties x=1; sin(x) ties π/3; exp(x) ties 0. Tikslios išvestinės yra [1,6], [0.5,−sqrt(3)/2], [1,1]. Lentelėje pateikti suapvalinti duomenys naudojami nepakeisti. Trečioji lentelė atitinka exp(−x).

`lagran()` interpoliuoja vidurinius 3 arba visus 5 taškus; `polyder()` randa daugianarių išvestines, `polyval()` jas įvertina ties viduriniu tašku.

## Integravimas

Pagal variantą integruojama x³−2x intervale [0,2], sin(x) intervale [0,π/2] arba exp(−x) intervale [0,1]. Tikslūs rezultatai: 0, 1, 1−exp(−1).

Rombergo metodui naudojamos log2(N)+1 eilutės, kad paskutinis trapecijų tinklelis turėtų N dalių.

Netiesioginio sin(x)/x integralo simbolinė reikšmė yra π/2. Metodai skaičiuoja intervale [0.001,100]. Todėl paklaida π/2 atžvilgiu apima ir praleistą pradžią bei uodegą. Atskira lentelės skiltis rodo tik skaitinę paklaidą baigtiniame intervale, kurio etaloną apskaičiuoja integral() su 1e-12 tolerancijomis. 20 taškų Gauss-Legendre tinklelis gali netiksliai atkurti ilgą svyruojantį intervalą; naudojamas būtent užduoties nurodytas taškų skaičius.
