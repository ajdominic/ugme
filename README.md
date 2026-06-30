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

**Citations**

If you use this repository or adapt the methods in your own work, please cite:

@article{dominic2023ugme,
  title = {Building insightful, memory-enriched models to capture long-time biochemical processes from short-time simulations},
  author = {Dominic, Anthony J. and Sayer, Thomas and Cao, Shixian and Markland, Thomas E. and Huang, Xuhui and Montoya-Castillo, Andrés},
  journal = {Proceedings of the National Academy of Sciences},
  volume = {120},
  number = {12},
  pages = {e2221048120},
  year = {2023}
}

and 

@article{dominic2023perspective,
  title = {Memory unlocks the future of biomolecular dynamics: Transformative tools to uncover physical insights accurately and efficiently}, 
  author = {Dominic, Anthony J. and Cao, Shixian and Huang, Xuhui and Montoya-Castillo, Andrés},
  journal = {Journal of the American Chemical Society},
  volume = {145},
  number = {18},
  pages = {9916−9927},
  year = {2023},
}
