# load_planetS.py
import requests
import os
from datetime import datetime
import xml.etree.ElementTree as ET
# import argparse
# import seaborn as sns
import pandas as pd


def check_internet_connection():
    """Check if there is internet access by pinging a known URL."""
    try:
        requests.get("https://www.google.com", timeout=5)
        return True
    except requests.ConnectionError:
        return False

def update_planets_exoplanet_catalog(catalog_url, output_file):
    """
    Updates the PlanetS exoplanet catalog from the DACE website.
    Warning: Work in progress, currently not working
    """
    # if True:
    #     print("PlanetS auto update currently not implemented, skipping updating attempt.")
    #     return

    if not check_internet_connection():
        print("No internet connection. Unable to update the exoplanet catalog.")
        return

    try:
        # Fetching the PlanetS catalog from internet. Code provided by Léna Parc + adapted by Artem Aguichine
        response = requests.get(catalog_url, verify=False, timeout=10)
        response.raise_for_status()  # Raise an error for HTTP issues

        # Parse XML content
        root = ET.fromstring(response.content)

        # VOTable namespace
        # This line defines a dictionary that maps the prefix 'v' so that the XML parser knows it belongs to the VOTable namespace defined by the IVOA standard
        ns = {'v': 'http://www.ivoa.net/xml/VOTable/v1.3'} 

        # Get field names (column headers)
        fields = root.findall(".//v:FIELD", ns)
        column_names = [f.attrib['name'] for f in fields]
        # print("Columns:", column_names)

        # Get data rows (in TABLEDATA format)
        rows = root.findall(".//v:TR", ns)

        # Extract data row-by-row
        data = []
        for tr in rows:
            values = [td.text for td in tr.findall("v:TD", ns)]
            data.append(values)

        # Convert to pandas dataframe
        df_planets_catalog_full = pd.DataFrame(data, columns=column_names)
        
        MJUP_TO_MEARTH = 317.828
        RJUP_TO_REARTH = 11.209
        df_planets_catalog_keep = pd.DataFrame({
            "pl_name": df_planets_catalog_full["Planet Name"],
            "pl_rade": pd.to_numeric(df_planets_catalog_full["Planet Radius [Rjup]"]) * RJUP_TO_REARTH,
            "pl_radeerr1": pd.to_numeric(df_planets_catalog_full["Planet Radius - Upper Unc [Rjup]"]) * RJUP_TO_REARTH,
            "pl_radeerr2": pd.to_numeric(df_planets_catalog_full["Planet Radius - Lower Unc [Rjup]"]) * RJUP_TO_REARTH,
            "pl_masse": pd.to_numeric(df_planets_catalog_full["Planet Mass [Mjup]"]) * MJUP_TO_MEARTH,
            "pl_masseerr1": pd.to_numeric(df_planets_catalog_full["Planet Mass - Upper Unc [Mjup]"]) * MJUP_TO_MEARTH,
            "pl_masseerr2": pd.to_numeric(df_planets_catalog_full["Planet Mass - Lower Unc [Mjup]"]) * MJUP_TO_MEARTH,
            "pl_eqt": pd.to_numeric(df_planets_catalog_full["Equilibrium Temperature [K]"]),
            "pl_peri": pd.to_numeric(df_planets_catalog_full["Orbital Period [days]"]),
            "pl_perierr1": pd.to_numeric(df_planets_catalog_full["Orbital Period - Upper Unc [days]"]),
            "pl_perierr2": pd.to_numeric(df_planets_catalog_full["Orbital Period - Lower Unc [days]"]),
            })

        # Writing the catalog to local file
        # Prepare the header
        header = [
            "# DACE PlanetS Catalog",
            f"# Source: {catalog_url}",
            f"# Catalog last updated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
        ]
        
        # Write the header and data to the file
        with open(output_file, "w") as f:
            f.write("\n".join(header) + "\n")

            f.write("# " + " ".join(df_planets_catalog_keep.columns) + "\n")

            df_planets_catalog_keep.to_csv(f, sep="\t", index=False, header=False)

        print(f"Catalog updated successfully and saved to {output_file}")
    except requests.RequestException as e:
        print(f"Error fetching the catalog: {e}")

def read_planets_last_update(output_file):
    """
    Reads the date of the last update from the catalog file and prints it.
    Parameters:
        output_file (str): Path to the catalog file.
    """
    if not os.path.exists(output_file):
        print("Catalog file does not exist.")
        return

    try:
        with open(output_file, "r") as f:
            for line in f:
                text_lookup = "# Catalog last updated:"
                if line.startswith(text_lookup):
                    print("# PlanetS catalog last updated: ",line.replace(text_lookup,'').replace('\n',''))
                    break
    except Exception as e:
        print(f"Error reading the catalog: {e}")

# File paths and URL of the PlanetS Catalog
planets_catalog_url = ("https://dace.unige.ch/downloads/exoplanets/PLANETS.vot")


output_file = "Models/catalog_exoplanets_planets.dat"
update_planets_exoplanet_catalog(planets_catalog_url, output_file)