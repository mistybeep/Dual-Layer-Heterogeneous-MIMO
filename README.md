# DBS-MA
This simulation code package is mainly used to reproduce the results of the following paper [1]:

[1] Yifan Guo, “Dual-end Fluid Antennas For Robust Anti-jamming in Low-altitude Air-ground Communications”

*********************************************************************************************************************************
If you use this simulation code package in any way, please cite the original paper [1] above. 
 
The author in charge of this simulation code pacakge is: Yifan Guo (email: guoyifan@nudt.edu.cn).

Please note that the MATLAB R2023a is used for this simulation code package,  and there may be some imcompatibility problems among different MATLAB versions. 

College of Electronic Science and Technology, National University of Defense Technology, Changsha 410073, China. 

*********************************************************************************************************************************
Abstract of the paper: 

This paper addresses the challenge of co-channel interference and intentional jamming in low-altitude air-ground communications.
Since conventional fixed-position antenna (FPA) systems lack spatial adaptability to dynamically balance signal enhancement against interference suppression, we propose a transformative fluid antenna system (FAS)-assisted heterogeneous dual-layer transmission architecture.
Specifically, a terrestrial base station with FPA serves ground users, while a low altitude-serving base station equipped with FAS communicates with the aerial user, also equipped with FAS, under the attack of a malicious jammer. 
We formulate a worst-case achievable rate maximization problem for aerial user subject to constraints including quality-of-service for terrestrial users, imperfect jamming directions, minimum antenna separation, etc.
To address the non-convex problem, we propose a fractional programming-block coordinate descent algorithm that alternately optimizes the transmit precoders, receive combiners, and antenna positions at both transceiver sides.
Convex hull-based approach and geometric boundary method are used to handle the jamming uncertainty and antenna placement constraints in confined spatial regions, respectively.
Extensive simulations validate significant performance gains. The FAS achieves up to 56\% higher data rates than FPA under equivalent power constraints. Strategic antenna repositioning demonstrably enhances signal quality while suppressing interference, maintaining robustness across diverse jammer channel uncertainties.

*********************************************************************************************************************************
Enjoy the reproducible research!