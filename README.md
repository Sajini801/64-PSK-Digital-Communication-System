# 64-PSK Digital Communication System 

## Overview
This project implements and analyses a 64 Phase Shift Keying (64 PSK) digital communication system using MATLAB.  
The system models the complete transmission chain, including vector encoding, modulation, transmission over an AWGN channel, demodulation and optimal detection at the receiver.
The project focuses on constellation based analysis to visualise the impact of noise and attenuation on signal integrity.


## Objectives
- Implement a complete 64 PSK modulation and demodulation chain
- Apply vector encoding and Gray code mapping to binary input data
- Model an AWGN channel with attenuation
- Perform Maximum Likelihood (ML) detection
- visualise transmitted and received constellation diagrams
- Analyse changes in symbol radius and noise effects

---

## Theory Background
64 PSK is a digital modulation technique where information is encoded in the **phase** of a carrier signal.  
Each symbol represents 6 bits of information. In an ideal channel, all symbols lie on a circle of fixed radius in the complex plane.  
In practical channels, attenuation reduces signal amplitude and noise introduces random changes, affecting detection performance.

---

## System Model
The communication system is modelled as:

- Transmitted symbol: s
- Channel attenuation: α
- Additive noise: n
- Received symbol: r = αs + n

This model displays amplitude reduction and noise corruption commonly observed in real world communication channels.

---

## Implemented Features
- Binary input recieving
- Vector encoding into 6bit symbols
- Gray code mapping
- 64 PSK modulation using complex exponentials
- AWGN channel simulation with attenuation
- Constellation plotting with consistent axis scaling
- Maximum Likelihood (ML) symbol detection
- Symbol and phase recovery

---

## Results and Observations
- The transmitted constellation forms a uniform circular pattern
- The received constellation exhibits a reduced radius due to attenuation
- Noise causes symbol spreading, increasing overlap between decision regions
- Using identical plot scaling highlights physical channel effects clearly
- ML detection provides optimal symbol decisions under AWGN conditions

---

## System Processing Stages

### Vector Encoding
The input binary stream is seperated into 6 bit blocks, with each block representing one symbol in the 64 PSK constellation.  
A predefined Gray code lookup table is used to map each 6 bit pattern to a symbol index from 0 to 63 which corresponds to a specific constellation point.  
For each 6 bit block, the matching Gray code is identified and the associated symbol index is assigned to the array, which is then used for modulation.

### Gray Coding
The symbol groups are Gray coded so that adjacent constellation points differ by only one bit, reducing the impact of symbol decision errors.

### Modulation (64 PSK)
Each Gray coded symbol index is mapped to a complex exponential with fixed magnitude and unique phase, producing a complex valued symbol array.

### Constellation Representation
The modulated symbols are represented in the I,Q plane, allowing visual inspection of phase separation and noise sensitivity.

### AWGN Channel
Complex Gaussian noise is added to each transmitted symbol, spreading the constellation points and modelling real world channel impairments.

### Demodulation
The received complex symbols are processed to estimate their phases and map them back to symbol indices.

### Maximum Likelihood (ML) Detection
Each received symbol is compared to all possible constellation points, and the symbol with the minimum Euclidean distance is selected.




