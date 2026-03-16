# Wireless-Comm.-of-RIS-assisted-Differential-Spatial-Modulation

[cite_start]This repository contains the MATLAB simulation code and the original research paper for the RIS-assisted Differential Transmitted Spatial Modulation (DTSM) scheme.

## 📖 Overview
The proposed DTSM scheme integrates the encoding process of differential spatial modulation (DSM) into a reconfigurable intelligent surface (RIS) assisted wireless communication system. At each time slot, only a single transmit antenna is activated to transmit an M-ary phase shift keying (PSK) symbol through the RIS.

**Key Advantages of DTSM:**
* **No CSI Requirement at the Receiver:** Due to the detection characteristics of DSM, the bit error rate (BER) performance remains satisfactory without requiring channel state information (CSI) estimation, thereby enhancing system robustness.
* **Shadow Area Mitigation:** The RIS adjusts the phase of the reflected signals, effectively mitigating shadow area fading and improving the signal-to-noise ratio (SNR) at the receiveR.
* **High Reliability in Complex Fading:** The scheme demonstrates robust BER performance across various scenarios, including Nakagami-$m$ fading channels.
* **Synchronization Independence:** Unlike Reflected Differential Spatial Modulation (RDSM), the proposed scheme eliminates the necessity for synchronization between the RIS and the receiver, while still reaping the benefits of DSM.

## 📝 Citation
If you find this repository helpful for your academic research, please consider citing our paper:

```bibtex
@article{zhang2024ris,
  title={RIS-assisted differential transmitted spatial modulation design},
  author={Zhang, Chaorong and Liu, Yuyan and Ng, Benjamin K. and Lam, Chan-Tong},
  journal={Signal Processing},
  year={2025},
  publisher={Elsevier},
  doi={10.1016/j.sigpro.2024.109767}
}

'''IEEE format:
C. Zhang, Y. Liu, B. K. Ng, and C.-T. Lam, “RIS-assisted differential transmitted spatial modulation design,” \textsl{Signal Process.}, vol. 230, p.109767, 2025.
