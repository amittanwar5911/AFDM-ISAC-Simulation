
# Automotive Radar Adaptation of AFDM-ISAC

## 1. Objective

The objective is to extend the existing AFDM-based Integrated Sensing and Communication (ISAC) simulation towards automotive radar applications.

The focus is on studying how AFDM-based sensing can support the detection and estimation of moving vehicles using range, Doppler, and target-detection processing.

The intended application is forward-looking automotive sensing in highway environments.

## 2. Existing Simulation Baseline

The current MATLAB implementation uses the following parameters:

| Parameter | Current Value |
|---|---:|
| Carrier frequency | 24 GHz |
| Number of subcarriers | 512 |
| Number of AFDM symbols | 32 |
| Subcarrier spacing | 22.729 kHz |
| Sampling frequency | 11.64 MHz |
| Chirp-periodic prefix | 64 samples |
| Radar targets | 3 |
| Modulation | 16-QAM |

### Existing Target Scenario

| Target | Range | Radial Velocity |
|---|---:|---:|
| Target 1 | 300 m | +20 m/s |
| Target 2 | 200 m | +8 m/s |
| Target 3 | 120 m | -8 m/s |

The simulation already implements:

- AFDM modulation and demodulation
- Delay-Doppler channel modelling
- Radar echo generation
- FCCR-based range-Doppler processing
- DAFT-domain sensing
- CA-CFAR target detection
- Range and velocity estimation

## 3. Automotive Radar Relevance

Automotive radar systems are used for applications such as:

- Forward collision warning
- Adaptive cruise control
- Vehicle detection and tracking
- Relative velocity estimation
- Blind-spot monitoring
- Advanced driver-assistance systems (ADAS)

The 76-81 GHz frequency range is particularly relevant to modern automotive radar applications.

## 4. Existing Limitations

The current simulation is a research prototype rather than an automotive-qualified radar implementation.

Its important limitations include:

1. The carrier frequency is currently 24 GHz rather than a dedicated 77-81 GHz automotive configuration.
2. The available bandwidth is approximately 11.64 MHz, resulting in a theoretical range resolution of approximately 12.89 metres.
3. The target scenario is not yet based on a moving ego vehicle and surrounding traffic.
4. Vehicle tracking, angle estimation and realistic road clutter are not implemented.
5. Automotive radar propagation and multipath conditions require further modelling.

## 5. Proposed Automotive Extensions

### Phase 1: Automotive Radar Parameter Configuration

- Introduce a configurable 77 GHz carrier-frequency scenario.
- Recalculate Doppler frequency and velocity ambiguity.
- Reassess the sampling frequency and signal bandwidth.
- Verify the range and velocity resolution.

### Phase 2: Highway Scenario Modelling

Introduce an automotive environment with:

- Ego vehicle
- Lead vehicle
- Adjacent-lane vehicle
- Oncoming or crossing vehicle
- Stationary roadside objects

### Phase 3: Improved Radar Sensing

Explore:

- Multi-target range and velocity estimation
- Improved CFAR detection
- Radar clutter and multipath modelling
- False-alarm and missed-detection analysis

### Phase 4: Automotive Sensing Evaluation

Evaluate:

- Range estimation accuracy
- Relative velocity estimation accuracy
- Detection probability
- False-alarm probability
- Performance under different SNR and Doppler conditions

## 6. Intended Outcome

The intended outcome is an automotive-oriented AFDM-ISAC simulation that investigates the potential of AFDM waveforms for joint wireless communication and vehicle radar sensing.

The project will provide a foundation for studying automotive radar signal processing, moving-target detection and high-mobility communication.

## 7. Future Research Direction

Further extensions may include:

- Automotive radar angle estimation
- MIMO radar modelling
- Vehicle trajectory tracking
- Radar-camera sensor fusion
- Comparison with conventional automotive radar waveforms

---

**Author:** Amit Tanwar  
**Project:** AFDM-Based Integrated Sensing and Communication
