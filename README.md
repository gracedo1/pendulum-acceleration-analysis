# Pendulum Acceleration Analysis Using Phyphox

Repository for the smartphone pendulum project conducted for General Physics. The experiment used Phyphox acceleration data to estimate the oscillation period of a pendulum made from a smartphone suspended by two parallel strings.

## Repository contents

```text
data/raw_xls/      Original Phyphox spreadsheet files
data/csv/          CSV files used for analysis
figures/           Figures generated from the analysis
results/           Trial-level and condition-level results
analysis.py        Python script used for period estimation
requirements.txt   Python package list
notebooks/         Colab/Jupyter notebook
```

## Experimental parameters

| Quantity | Value |
|---|---:|
| Pendulum length, L | 1.20 m |
| Upper support separation | 18 cm |
| Envelope hole separation | 18 cm |
| Smartphone + envelope mass | 239 g |
| Smartphone + envelope + added mass | 445 g |
| Smartphone height increase with cardboard | 8.5 cm |
| Gravitational acceleration used | 9.80 m/s^2 |

The small-angle reference period was calculated as

```text
T = 2*pi*sqrt(L/g) = 2.199 s
```

## Experimental conditions

| File label | Condition |
|---|---|
| angle10_trial1-3 | Basic pendulum, 10 degree initial angle |
| angle20_trial1-3 | Basic pendulum, 20 degree initial angle |
| angle30_trial1-3 | Basic pendulum, 30 degree initial angle |
| cardboard20_trial1-3 | 20 degree initial angle, smartphone raised with cardboard |
| mass20_trial1-3 | 20 degree initial angle, added mass |

## Running the analysis

Install the required packages and run the analysis script from the repository root.

```bash
pip install -r requirements.txt
python analysis.py
```

The script reads the CSV data in `data/csv/`, estimates periods from the y-axis acceleration signal, and writes the output tables and figures to `results/` and `figures/`.

## Summary of results

| Condition | Mean period (s) | Standard deviation (s) | Mean error (%) |
|---|---:|---:|---:|
| 10 deg | 2.225 | 0.075 | 2.873 |
| 20 deg | 2.262 | 0.010 | 2.863 |
| 30 deg | 2.263 | 0.025 | 2.918 |
| 20 deg + cardboard | 2.058 | 0.028 | 6.398 |
| 20 deg + added mass | 2.254 | 0.035 | 2.506 |

The basic 10, 20, and 30 degree trials gave similar periods. The added-mass trials remained close to the basic 20 degree condition, while the cardboard condition gave a shorter period because the smartphone position changed the effective pendulum length.

## Reference

R. Mathevet, N. Lamrani, L. Martin, P. Ferrand, J. P. Castro, P. Marchou, and C. M. Fabre, “Quantitative analysis of a smartphone pendulum beyond linear approximation: A lockdown practical homework,” *American Journal of Physics*, vol. 90, no. 5, pp. 344–350, 2022.
