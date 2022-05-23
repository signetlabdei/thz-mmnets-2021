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
- Compatibility with MPC channel models


Model Description, Usage and Validation of the TeraSim
------------------------------------------------
Please refer to the thz.rst in the thz/doc folder



MPC Channel model usage Example
------------------------------------------------
Terasim supports two distinct kinds of channel models [1][2]. check [3] for more information regarding implementation.
macroMpcChannel is an illustration of MPC channel model usage.
[2] The channel model referred to as "FS" in the code is fully stochastic and can be utilized easily.
[1] refer as "HB" in code is part of Raytracing, and these steps must be taken prior to using the module:

download and install "qd-channel" module(https://github.com/AmirAshtariG/qd-channel.git).
to generate Raytracing:
1) current examples can be found in "qd-channel/model/QD/".
or
2) "qd-realization" (https://github.com/signetlabdei/qd-realization) or any other Raytracing tool ("qd-realization" has already ns-3 compatible output for "qd-channel" module) can be used to generate new scenarios.

[1]Yi Chen, Yuanbo Li, Chong Han, Ziming Yu, and Guangjian Wang. 2021. Channel
Measurement and Ray-Tracing-Statistical Hybrid Modeling for Low-Terahertz Indoor Communications. IEEE Transactions on Wireless Communications ,Early Access (2021). https://doi.org/10.1109/TWC.2021.3090781
[2]Shihao Ju, Yunchou Xing, Ojas Kanhere, and Theodore S. Rappaport. 2021. Millimeter Wave and Sub-Terahertz Spatial Statistical Channel Model for an Indoor Office Building. IEEE Journal on Selected Areas in Communications 39, 6 (Apr2021), 1561–1575.
[3]A. Ashtari Gargari, M. Polese, M. Zorzi, "Full-stack comparison of channel models for networks above 100 GHz in an indoor scenario", in Proceedings of the 5th ACM Workshop on Millimeter-Wave and Terahertz Networks and Sensing Systems (mmNets '21), October 2021, Pages 43–48, https://doi.org/10.1145/3477081.3481677


Limitation
-----------------------------------
Due to support for MPC channel models and Terasim's reliance on the qd-channel module, Terasim requires installation of the qd-channel module despite not employing it.

Copy Right
------------------------------------------------
Northeastern University https://unlab.tech/
