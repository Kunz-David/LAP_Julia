[docs-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-dev-url]: https://kunz-david.github.io/LAP_Julia.jl/dev/
[documentation-build-img]: https://github.com/Kunz-David/LAP_Julia.jl/workflows/Documentation/badge.svg

# LAP_Julia

_Local all-pass filtering registration method implemented in Julia from [this paper](http://www.ee.cuhk.edu.hk/~tblu/monsite/pdfs/gilliam1701.pdf)._

| **Documentation**                                  | **Build Status**                                  |
|:--------------------------------------------------:|:-------------------------------------------------:|
| [![][docs-dev-img]][docs-dev-url]                  |  ![Documentation][documentation-build-img]        |


_Check out the new [documentation][docs-dev-url]!_



## Installation

For the latest version of this module open up a Julia terminal and type:
```Julia
using Pkg; Pkg.add(PackageSpec(url="https://github.com/Kunz-David/LAP_julia"))
using LAP_julia
```

## Stav práce


### 22.7.2020
- Podarilo se mi docilit dalsich zrychleni na implementaci LAP a PF-LAP lap metod, bohuzel je `imfilter` funkce na filtraci
obrazku v Julii pomalejsi nez `imfilter` v Matlabu, s tim nic moc neudelam, muzu predvest v exprimentu.
- Zrychlil jsem implementaci obou sparse metod: `sparse_pflap` a `sparse_lap` a pridal jsem _fitovani na globalni kvadraticky model_.
  - Obe vyber fitovani a dalsich paramatru spusteni je v dokumentaci.
  - Fitovani na kvadraticky model je rychle a pro konstantni displacement pole je presnejsi nez PF-LAP viz priklad.
- pridano nekolik testovacich funkci `test_registration_alg`, `asses_flow_quality` ...
- nekolik experimentu: viz `test/experiments/...`
  - mezi residualem linearnich rovnic LAP metody a chybou odhadu posunu je korelace - prumerna `Normalized Cross-correlation(ncc)` = `0.22` (v experimentu je histogram rozlozeni) - nevim jak toho vyuzit zatim, musi se to docela narocne zpetne spocitat.
- snazim se zprovoznit benchmark BIRL pro moje methody.


__Test rychlosti:__
- obrazek: `257x257`, displacement pole: konstatni, max posun: `15 pix`, (`img, imgw, flow = gen_init(:lena, :uniform, flow_args=[-1 - 4im, 15]);`), pocet mereni: `10`

__sparse_pflap__:
- ARG metody: pocet bodu: `459`, max povolena vzdalenost bodu: `10 pix`
```julia
timer=TimerOutput("sparse pf lap");
method_kwargs = Dict(:timer => timer, :display => false, :max_repeats => 1, :point_count => 500, :spacing => 10)
flow_est, source_reg, timer, results = test_registration_alg(sparse_pflap, img, imgw, flow, [], method_kwargs, timer=timer)
```

_sparse_pflap_ - rychlost:
```
────────────────────────────────────────────────────────────────────────────────────────
                                                Time                   Allocations      
                                        ──────────────────────   ───────────────────────
           Tot / % measured:                 34.0s / 16.4%           3.98GiB / 98.2%    

Section                         ncalls     time   %tot     avg     alloc   %tot      avg
────────────────────────────────────────────────────────────────────────────────────────
sparse pf lap                       10    5.57s   100%   557ms   3.91GiB  100%    401MiB
  single filter pyramid level       70    5.17s  92.7%  73.8ms   3.85GiB  98.3%  56.3MiB
    single lap at points            70    4.12s  73.9%  58.9ms   3.32GiB  84.7%  48.5MiB
      filtering                     70    3.00s  53.9%  42.9ms    687MiB  17.1%  9.81MiB
      prepare A and b               70    551ms  9.88%  7.87ms    417MiB  10.4%  5.96MiB
        window sum part 1          210    303ms  5.43%  1.44ms    250MiB  6.25%  1.19MiB
        window sum part 2          140    247ms  4.43%  1.76ms    166MiB  4.15%  1.19MiB
      multi mat div                 70    543ms  9.74%  7.75ms   2.13GiB  54.3%  31.1MiB
      calculate flow                70    545μs  0.01%  7.78μs   1.01MiB  0.03%  14.8KiB
    interpolate image               70    824ms  14.8%  11.8ms    420MiB  10.5%  6.00MiB
    interpolate flow                70    189ms  3.38%  2.69ms   88.8MiB  2.22%  1.27MiB
  setup                             10    407ms  7.31%  40.7ms   67.1MiB  1.67%  6.71MiB
    find edge points                10    285ms  5.11%  28.5ms   41.9MiB  1.04%  4.19MiB
    hist match                      10   88.4ms  1.59%  8.84ms   5.15MiB  0.13%   528KiB
────────────────────────────────────────────────────────────────────────────────────────
```

_sparse_pflap_ - presnost:
```
mse            | 0.001
rmse           | 0.037
time           | 10.343
ncc            | 0.981
flow_mae       | 0.243
angle-rmse     | 0.592
angle-mae      | 0.489
mae            | 0.019
flow_rmse      | 0.09
```

__pflap__:

```julia
timer=TimerOutput("pf lap");
method_kwargs =Dict(:timer => timer, :display => false, :max_repeats => 1)
flow_est, source_reg, timer, results = test_registration_alg(polyfilter_lap, img, imgw, flow, [], method_kwargs, timer=timer);
```

_pflap_ - rychlost:
```
────────────────────────────────────────────────────────────────────────────────────────
                                                Time                   Allocations      
                                        ──────────────────────   ───────────────────────
           Tot / % measured:                 16.9s / 63.3%           4.80GiB / 98.7%    

Section                         ncalls     time   %tot     avg     alloc   %tot      avg
────────────────────────────────────────────────────────────────────────────────────────
pf lap                              10    10.7s   100%   1.07s   4.73GiB  100%    484MiB
  single filter pyramid level       70    10.3s  95.9%   147ms   4.68GiB  99.0%  68.5MiB
    LAP                             70    5.97s  55.6%  85.2ms   2.85GiB  60.3%  41.7MiB
      filtering                     70    3.21s  29.9%  45.8ms    687MiB  14.2%  9.81MiB
      multli mat div gem            70    1.04s  9.68%  14.8ms    840MiB  17.3%  12.0MiB
      prepare A and b               70    957ms  8.93%  13.7ms    521MiB  10.8%  7.45MiB
        window sum part 1          210    453ms  4.23%  2.16ms    250MiB  5.16%  1.19MiB
        window sum part 2          140    301ms  2.81%  2.15ms    166MiB  3.44%  1.19MiB
      calculate flow                70    152ms  1.42%  2.18ms    140MiB  2.89%  2.00MiB
    inpainting                      70    1.79s  16.7%  25.6ms   1.10GiB  23.2%  16.1MiB
      replicating borders           70   82.3ms  0.77%  1.18ms    124MiB  2.56%  1.77MiB
    smoothing                       70    1.44s  13.4%  20.5ms    258MiB  5.32%  3.68MiB
    image interpolation             70    989ms  9.22%  14.1ms    420MiB  8.67%  6.00MiB
  setup                             10    439ms  4.09%  43.9ms   50.0MiB  1.03%  5.00MiB
────────────────────────────────────────────────────────────────────────────────────────
```

_pflap_ - presnost:
```
mse            | 0.003
rmse           | 0.052
time           | 10.725
ncc            | 0.989
flow_mae       | 0.54
angle-rmse     | 3.13
angle-mae      | 1.033
mae            | 0.04
flow_rmse      | 2.126
```

### 9.7.2020
- Julia implementace PF-LAP a LAP metod je ted o neco rychlejsi nez implementace v Matlabu
  - PF-LAP Matlab implementace cca 1.5s na obrazky 256x256, Julia cca 1s
  - LAP Matlab implementace cca 0.1s na obrazky 256x256, Julia cca 0.09s
  - __TODO: zrychleni sparse metod__

__Speedtest:__
```
────────────────────────────────────────────────────────────────────────────────────────
                                                Time                   Allocations      
                                        ──────────────────────   ───────────────────────
           Tot / % measured:                 1.05s / 100%             549MiB / 100%     

Section                         ncalls     time   %tot     avg     alloc   %tot      avg
────────────────────────────────────────────────────────────────────────────────────────
polyfilter lap                       1    1.04s   100%   1.04s    549MiB  100%    549MiB
  single filter pyramid level        7    1.04s   100%   149ms    545MiB  99.4%  77.9MiB
    single lap                       7    642ms  61.5%  91.7ms    332MiB  60.5%  47.4MiB
      filtering                      7    293ms  28.0%  41.8ms   35.9MiB  6.54%  5.12MiB
      prepare A and b                7    139ms  13.3%  19.8ms   75.4MiB  13.7%  10.8MiB
        window sum 1                21   82.6ms  7.91%  3.93ms   39.0MiB  7.10%  1.85MiB
        window sum 2                14   26.7ms  2.56%  1.91ms   26.0MiB  4.73%  1.85MiB
      multli mat div                 7   78.0ms  7.47%  11.1ms   84.0MiB  15.3%  12.0MiB
      calculate flow                 7   61.4ms  5.88%  8.77ms   49.0MiB  8.94%  7.01MiB
    inpainting                       7    202ms  19.4%  28.9ms    135MiB  24.6%  19.3MiB
    smoothing                        7    122ms  11.7%  17.5ms   32.7MiB  5.97%  4.68MiB
    image interpolation              7   69.6ms  6.67%  9.94ms   38.5MiB  7.02%  5.50MiB
  setup                              1   2.91ms  0.28%  2.91ms   3.50MiB  0.64%  3.50MiB
────────────────────────────────────────────────────────────────────────────────────────
```
- Pridana moznost casovat registracni metody.
_Priklad:_

```julia
using TimerOutputs, LAP_julia

TimerOutputs.enable_debug_timings(LAP_julia)
timer = TimerOutput("Registration");
@timeit timer "polyfilter lap" begin
    flow_est, source_reg = polyfilter_lap(img, imgw, display=false, timer=timer)
end
print_timer(timer)
```
- Novy zpusob generovani pixeloveho posunu.
  - pripraveno na testovani fitovani na kvadraticky polynom
  - funkce `gen_quad_flow`
- TODO: kvadraticka interpolace posunu (fitovani).


### 29.5.2020
- projektova zprava hotova
- porovanani s metodami v Matlabu nedopadly dobre, moje Julia implementace je cca 5x pomalejsi
  - TODO: profilovani metod, zjistit co se deje
- TODO: dokumentace sparse metod
- pridat doklikatelne Jupyter experimenty.


### Shrnuti stavu ke dni 18.5.2020

#### Stav
- Jak testuju funkcnost meho kodu:
  - Vytvoreni testovanych obrazku; original a deformovany obrazek
    - Generovani nahodne hladke deformace pro testovani je implementovano v `gen_tiled_flow`
      - tady bych chtel pridat jeste jiny zpusob vytvoreni deformace, protoze deformace vypada hezky :) (viz dokumentace), ale nevim jak vypada deformace opravdovych obrazku. Napriklad [tady](https://www.researchgate.net/publication/316876842_Iterative_fitting_after_elastic_registration_An_efficient_strategy_for_accurate_estimation_of_parametric_deformations) je defomace jednodussi a myslim, ze by nase metoda mohla fungovat na jednodussi a pomalejsi deformace fungovat lepe. Chtel bych na to pouzit [tuto](https://evizero.github.io/Augmentor.jl/) julia knihovnu.
    - Deformovany obrazek ziskam z originalu za pouziti funkce `warp_img`, ta interpoluje linearne. (Myslim, ze se nevyplati interpolovat jinak, protoze presnost se tolik nezmeni, ale rychlost je znatelne mensi. ??)

- Metody registrace:
  - Funkce `single_lap` z paperu funguje dle ocekavani. (viz dokumentace)
  - Funkce `polyfilter_lap` z paperu funguje dle ocekavani. (viz dokumentace)
  - Funkce `single_lap_at_points`
    - funguje dobre -> na obrazek 256x256, 25 klicovych bodu, deformovan pomalu se menici deformaci a s malym maximalnim posunutim,
    je cca 6x rychlejsi nez single_lap a podobne presna, pokud jsou body dobre distribuovane.
    - spolecne s vyberem bodu a interpolaci s `interpolate_flow` je to cele 5x rychlejsi. Dole pridavam kod na vyzkouseni.
    - je potreba zjistit za jakych podminek to funguje nejlepe. Udelat nejake testy??
    - je potreba dodelat dokumentaci.
  - Funkce `sparse_pflap` je zatim stejne rychla jako `polyfilter_lap` a je o neco nepresnejsi

- Interpolace/fitovani globalni deformace:
  - je implementovano pomoci RBF interpolace ve funkci `interpolate_flow`.
  - mam v planu zkusit jednodussi globalni fit funkci, treba kvadraticky polynom.

- Dokumentace:
  - obsahuje dokumentace vsech vyznamnych public funkci a vetsinu privatnich, krom single_lap_at_points.
  - obsahuje navod na volani metod `polyfilter_lap` a `single_lap`, jejich porovnani a ukazky vystupu.
  - je potreba dodelat Examples do hezci formy.
  - chtel bych jeste pridat nejake lepsi porovnani vysledku metody a originalu, asi by stacilo neco jako `showflow(orig .- estim)`

#### Plan praci
1) Dodelat `sparse_pflap` aby fungovala presnejsi
2) Psat psat psat (potrebuji mit 10 normostran do 29.5.2020) -> formou na semestralni projekt
  a) co je registrace, na co je a porovnani s tim co se dela ted ve svete
  b) z jakych metod vychazim a jak se je snazim zrychlit/zlepsit
  c) jak to probihalo, co jsem pouzil, co jsem zkusil...
  d) jake jsou moje vysledky
3) male prace -> dokumentace?
4) ...


Kod pro vyzkouseni `single_lap_at_points`
```julia
using LAP_julia
include("useful.jl")

img = testimage("lena_gray")
img = Float32.(img)

flow = gen_tiled_flow(size(img), 15, 100)
showflow(flow)

# generate warpped image
imgw = LAP_julia.interpolation.warp_img(img, -real(flow), -imag(flow));

# params:
#lap
filter_num = 3;
fhs = 25;
window_size = [51, 51];
#points
point_count = 25;
spacing = 40;

mask = parent(padarray(trues(size(img).-(2*fhs, 2*fhs)), Fill(false, (fhs, fhs), (fhs, fhs))))

# get points
inds = find_edge_points(img, sigma = 1, number = point_count, spacing = spacing, mask = mask);

# reform points:
pos_x = [ind[1] for ind in inds]
pos_y = [ind[2] for ind in inds]
points = LAP_julia.inds_to_points(inds)

# see points
imgshow(img)
PyPlot.scatter(pos_y, pos_x, marker = :x); gcf()

# run methods
flow_est_all = single_lap(img, imgw, fhs, window_size, 3)
# flow_est_points = single_lap(img, imgw, fhs, window_size, 3, points)
flow_est_new, all_coeffs = single_lap_at_points(img, imgw, fhs, window_size, 3, points)

completed_flow = interpolate_flow(flow_est_new, inds)

# compare resulting estimations:
# show calculated flows
showflow(completed_flow, figtitle="flow estimated from lap at points")
showflow(flow, figtitle="original flow")

# show difference and mean squared error of both
showflow(flow .- completed_flow, figtitle="orig - estim from points")
LAP_julia.mse(flow, completed_flow)
mag_mse_points = LAP_julia.vec_len(LAP_julia.mse(flow, completed_flow))
println("Points: magnitude of the vector of the mean squared error is: ", mag_mse_points)

showflow(flow .- flow_est_all, figtitle="orig - estim classic single")
inpainted = copy(flow_est_all); LAP_julia.inpaint_nans!(inpainted)
LAP_julia.mse(flow, inpainted)
mag_mse_classic = LAP_julia.vec_len(LAP_julia.mse(flow, inpainted))
println("Classic: magnitude of the vector of the mean squared error is: ", mag_mse_classic)

# whole process speed comparison:
# new method
bench_find_points = @benchmark inds = find_edge_points(img, sigma = 1, number = point_count, spacing = spacing, mask = mask)
bench_lap_points = @benchmark flow_est_new, all_coeffs = single_lap_at_points(img, imgw, fhs, window_size, 3, points)
bench_fit_points = @benchmark completed_flow = interpolate_flow(flow_est_new, inds)

# classic single lap
bench_single_lap = @benchmark flow_est_all = single_lap(img, imgw, fhs, window_size, 3)

# speedup
new_speed = median(bench_find_points.times) + median(bench_lap_points.times) + median(bench_fit_points.times)
classic_speed = median(bench_single_lap.times)
speedup = classic_speed/new_speed

println("The speedup is: " speedup)
```
