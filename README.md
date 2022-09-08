# Tw-STI
Code given in support of ASIACRYPT 2022 paper, available on the IACR ePrint Archive [ePrint:2021/1384](https://eprint.iacr.org/2021/1384):

> **Log-S-unit Lattices Using Explicit Stickelberger Generators to Solve Approx Ideal-SVP**
> Olivier Bernard, Andrea Lesavourey, Tuong-Huy Nguyen and Adeline Roux-Langlois

This work is supported by the European Union PROMETHEUS project (Horizon 2020 Research and Innovation Program, grant 780701).


### Environment

The provided code has been tested with the following versions.

| Soft      | Version  | Url  |
| :---      | :---     | :--- |
| Magma     | v2.24-10 | https://magma.maths.usyd.edu.au/magma/ |
| SageMath  | v9.0     | https://www.sagemath.org/ |
| fplll     | v5.3.2   | https://github.com/fplll/fplll |


##### Remark
We suppose that fplll is in `/usr/local/bin`.
For changing this, you must edit in `./src/lattice.py` the following line:
```
__FPLLL_PATH = "/usr/local/bin/";
```

### Computational workflow

Note that `./data/list_ms_pcmp` contains the list of cyclotomic fields conductors for which h+=1 is known/computable, and of degree 21 < n < 211.

For simplicity, we suppose that everything is executed from `./scripts`, but it should work from everywhere. The conductor of the targetted cyclotomic field is `<m>`. Each script accepts a list of conductors, but beware that one or several threads will be launched for **each** of the conductors in the list.

0. Make sure to create a `logs` folder besides the `scripts` folder.
1. Compute places: complex embeddings, factor bases:
```
./cf_places.sh <m>
```
2. Compute circular units / real relative norm generators / Stickelberger generators:
```
./urs.sh <m>
```
3. Saturate (2 saturation, one pass):
```
./saturation.sh <m>
```
4. Compute S-units; for this Magma is needed, and it should work up to degree n < 80:
```
./sunits.sh <m>
```
5. Compute log-S-unit (sub)lattices associated to these family (urs/sat/[su]) for each of the {iso/noiso}x{exp/tw} options
```
./sti_twphs_lattice.sh <m>
```
6. Reduce each available lattice with LLL and BKZ-40
```
./lat_reduce.sh <m>
```
7. Precompute Gram-Schmidt orthogonalizations:
```
./gso.sh <m>
```
8. Evaluate the geometry of the obtained lattices: root hermite factor, orthogonality defect, table of Gram-Schmidt log norms.
```
./eval_geo.sh <m>
```
9. Simulate approximation factors on 100 random targets for split prime ideals of size 2^100, and also run CDW on these targets:
```
./rand_targets.sh <m>
./approx_factor.sh <m>
./cdw_protocol.sh <m>
```


### Plotting results

First, create a `./figures` folder besides the `./scripts` folder.

##### Obtaining Fig. 1.1, 5.3 and 5.4 
Run the `plot_preproc.sage` and then the `plot_minGF.py` script:
```
./plot_preproc.sage
./plot_minGF.py
```
This creates files `GF_CDW-mprandom2_intro.png` (Fig.1.1), `GF_CDW-mprandom2.png` (Fig.5.2) and `zoomGF_with_cdw_lower-mprandom2.png` (Fig.5.3) in folder `./figures`.


##### Obtaining Fig. 5.1, C.1-3
Run the `plot_raw_bkz_2fig.py` script with two chosen conductors m1 and m2 for orb orbits as:
```
./plot_raw_bkz_2fig.py <orb> <m1> <m2>
```
This creates file `z<m1>-z<m2>_comparison_raw_bkz_d<orb>.png` in folder `./figures`.


##### Obtaining Fig. C.4
Run the `plot_sets_2fig.py` script with two chosen conductors m1 and m2 for orb orbits as:
```
./plot_sets_2fig.py <orb> <m1> <m2>
```
This creates file `z<m1>-z<m2>_comparison_sets_d<orb>.png` in folder `./figures`.



### Files organisation

##### Folder list
- `./src`: All interesting stuff.
- `./scripts`: Sage/Magma scripts. Each script shall be called _via_ its Bash `.sh` counterpart.
This will redirect logs, detach thread(s), detect the number of orbits and the available S-unit sets.
- `./data`: this is where all data are put. Beware that the whole precomputation for the 192 fields weights > 1To.
- `./logs`: logs of computations, including timings and many diagnostics.
- `./figures`: plotted results


##### Naming conventions for precomputations


|Extension | Content |
|:---|:---|
`.inf`| Complex embeddings
`.fb`|  Factor bases (from d=1 to d=dmax predicted by Twisted-PHS)
`.targets`| Random targets, in the form of infinite absolute values, and p-adic valuations
`.urs`| Circular units, Real S+-units (h+=1), Stickelberger generators
`.sat`| 2-Saturated (one pass) family from urs 
`.su`| S-units, complete group when available (up to degree 80)
`.lat`,`.lat_gso`| Lattice diretly obtained after Twisted-PHS-like construction 
`.lll`, `.lll_U`,`.lll_gso`| LLL Version of the above, with unitary transformation
`.bkz-<bk>`, `.bkz-<bk>_U`,`.bkz-<bk>_gso`| BKZ reduction of block size bk, with unitary transformation and Gram-Schmidt Orthogonalization
`.geo`| Tables of all geometric indicators for this conductor
`.gsn`| Gram-Schmidt log norms (each gsn file contains data for raw/lll/bkz lattice for one given option)
`.afinf`,`.gf`,`.afsup`,`.hf`| Estimated approximation factors (Minkowski, Gaussian Heuristic, AGM inequality, Hermite Factor



### License

&copy; 2021, Olivier Bernard, Andrea Lesavourey, Tuong-Huy Nguyen and Adeline Roux-Langlois.  
This work is published under the GNU General Public License (GPL) v3.0.  
See the LICENSE file for complete statement.

