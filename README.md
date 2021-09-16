Channel Modeling for TeraSim
------------------------------------------------
This repo features the implementation of two additional channel models for TeraSim (see below). The implementation is based on the paper

> A. Ashtari Gargari, M. Polese, M. Zorzi, "Full-stack comparison of channel models 
> for networks above 100 GHz in an indoor scenario", in Proceedings of the 5th ACM Workshop 
> on Millimeter-Wave and Terahertz Networks and Sensing Systems (mmNets '21), October 2021, 
> Pages 43–48, https://doi.org/10.1145/3477081.3481677

Please cite the paper ([bibtex](https://dl.acm.org/doi/10.1145/3477081.3481677)) if you use this code.


TeraSim
------------------------------------------------
TeraSim is the first simulation platform for THz communication networks which captures the capabilities of THz devices and the peculiarities of the THz channel. The simulator has been developed considering two major types of application scenarios, namely, nanoscale communication networks (average transmission range usually below one meter) and macroscale communication networks (distances larger than one meter). The simulator consists of a common channel module, separate physical and link layers for each scenario, and two assisting modules, namely, THz antenna module and energy harvesting module, originally designed for the macroscale and nanoscale scenario, respectively. TeraSim is expected to enable the networking community to test THz networking protocols without having to delve into the channel and physical layers. 

A detailed explanation of the module can be found in [TeraSim: An ns-3 extension to simulate Terahertz-band communication networks](https://doi.org/10.1016/j.nancom.2018.08.001).


v1.1 release notes
------------------------------------------------
The latest TeraSim v1.1 release includes the following improvements:
- Implementation of [ADAPT: An Adaptive Directional Antenna Protocol for medium access control in Terahertz communication networks](https://doi.org/10.1016/j.adhoc.2021.102540). This is a MAC protocol for the macroscale scenario.
- Optimization of the code for improved computational efficiency
- Compatibility with ns-3.33


Model Description, Usage and Validation of the TeraSim
------------------------------------------------
Please refer to the thz.rst in the thz/doc folder


Copy Right
------------------------------------------------
Northeastern University https://unlab.tech/
