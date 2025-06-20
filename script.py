import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.ticker import FuncFormatter

# Load Data
point_loma = pd.read_csv('/Users/elizabethmurphy/Downloads/PointLoma_sewage_qPCR.csv', parse_dates=["Sample_Date"])
south_bay = pd.read_csv('/Users/elizabethmurphy/Downloads/SouthBay_sewage_qPCR.csv', parse_dates=["Sample_Date"])
encina = pd.read_csv('/Users/elizabethmurphy/Downloads/Encina_sewage_qPCR.csv', parse_dates=["Sample_Date"])

# Add sewershed labels
point_loma["Sewershed"] = "Point Loma"
south_bay["Sewershed"] = "South Bay"
encina["Sewershed"] = "Encina"

# Combine datasets
data = pd.concat([point_loma, south_bay, encina], ignore_index=True)

# 2023 (all sewersheds)
mask_2023 = (data["Sample_Date"] >= "2023-01-01") & (data["Sample_Date"] <= "2023-12-31")
data_2023 = data.loc[mask_2023]

sns.lineplot(data=data_2023, x="Sample_Date", y="Mean viral gene copies/L", hue="Sewershed", marker="o")
plt.title("SARS-CoV-2 Viral Load by Sewershed (2023)")
plt.xlabel("Sample Date")
plt.ylabel("Viral Load (gene copies/L)")
plt.legend(title="Sewershed")
plt.show()

# 2024 (all sewersheds)
mask_2024 = (data["Sample_Date"] >= "2024-01-01") & (data["Sample_Date"] <= "2024-12-31")
data_2024 = data.loc[mask_2024]

sns.lineplot(data=data_2024, x="Sample_Date", y="Mean viral gene copies/L", hue="Sewershed", marker="o")
plt.title("SARS-CoV-2 Viral Load by Sewershed (2024)")
plt.xlabel("Sample Date")
plt.ylabel("Viral Load (gene copies/L)")
plt.legend(title="Sewershed")
plt.show()

# June–September 2022 (all sewersheds)
mask = (data["Sample_Date"] >= "2022-06-01") & (data["Sample_Date"] <= "2022-09-30")
data_filtered = data.loc[mask]

sns.lineplot(data=data_filtered, x="Sample_Date", y="Mean viral gene copies/L", hue="Sewershed", marker="o")
plt.title("SARS-CoV-2 Viral Load by Sewershed (June–Sept 2022)")
plt.xlabel("Sample Date")
plt.ylabel("Viral Load (gene copies/L)")
plt.legend(title="Sewershed")
plt.show()
