## UGME

The `ugme` is a software for building and analyzing time-convolutionless generalized master equation models of biomolecular dynamics. We developed this method to directly include short-time memory effects in reduced dynamical models that recover long-time kinetics from short-time simulations. 

**Overview**

Our U-GME framework provides a way to incorporate non-Markovian effects into reduced dynamical models while retaining a compact description of the dynamics. We further introduce a simple averaging procedure to tame the noise from underconverged correlation matrices.

This repository contains tools for:
- Estimating time-local memory functions from a time-series of transition probability matrices (TPMs), $C(t)$
- Computing the time-local generator $U(t)$ from the $C(t)$
- Characterizing memory plateau $\tau_R$ timescales with the RMSE error metric (in the absence of noise)
- Identifying optimal onset ($t_r$) and offset $(\tau_R)$ of averaging for the time-local propagator U(t)
- Propagating the long-time dynamics of the TPMs with the resulting $U_\infty(\delta t) \equiv U(t \geq \tau_R)$.

## Associated Publications

The methods and code in this repository were developed in connection with the following publications:

1. Dominic, A. J. III; Sayer, T.; Cao, S.; Markland, T. E.; Huang, X.; Montoya-Castillo, A.  
   **Building insightful, memory-enriched models to capture long-time biochemical processes from short-time simulations.**  
   *Proceedings of the National Academy of Sciences* **2023**, *120* (12), e2221048120.

2. Dominic, A. J. III; Cao, S.; Huang, X.; Montoya-Castillo, A.  
   **Memory unlocks the future of biomolecular dynamics: Transformative tools to uncover physical insights accurately and efficiently.**  
   *Journal of the American Chemical Society* **2023**, *145* (18), 9916–9927.

## Citation

If you use this repository or adapt the methods in your own work, please cite:

```bibtex
@article{dominic2023ugme,
  title = {Building insightful, memory-enriched models to capture long-time biochemical processes from short-time simulations},
  author = {Anthony J. Dominic III and Thomas Sayer and Siqin Cao and Thomas E. Markland and Xuhui Huang and Andrés Montoya-Castillo},
  journal = {Proceedings of the National Academy of Sciences},
  volume = {120},
  number = {12},
  pages = {e2221048120},
  year = {2023}
}
```

and 

```bibtex
@article{dominic2023perspective,
  title = {Memory unlocks the future of biomolecular dynamics: Transformative tools to uncover physical insights accurately and efficiently}, 
  author = {Anthony J. Dominic III and Siqin Cao and Xuhui Huang and Andrés Montoya-Castillo},
  journal = {Journal of the American Chemical Society},
  volume = {145},
  number = {18},
  pages = {9916−-9927},
  year = {2023},
}
```
