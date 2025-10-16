import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path


df = pd.read_csv('6_14_10_mp03b.csv')

# ---- Plot ----
plt.figure(figsize=(7, 5))
plt.plot(df["Temperature"], df["kappa_ph"],   label=r"$\kappa_{\mathrm{ph}}$")
plt.plot(df["Temperature"], df["kappa_diff"], label=r"$\kappa_{\mathrm{diff}}$")
plt.plot(df["Temperature"], df["Kappa_total"],label=r"$\kappa_{\mathrm{total}}$")

plt.xlabel("Temperature (K)")
plt.ylabel(r"Thermal conductivity $\kappa$ (W m$^{-1}$ K$^{-1}$)")
plt.legend()
plt.tight_layout()

plt.savefig("kappa_vs_T.png")

plt.show()

