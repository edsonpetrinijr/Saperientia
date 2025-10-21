import csv
import numpy as np

def parse_ra_dec(ra_hrs, ra_min, ra_sec, dec_sign, dec_deg, dec_min, dec_sec):
    """Convert RA/Dec from sexagesimal to decimal degrees."""
    # RA: hours to degrees (15 degrees per hour)
    ra_decimal = (float(ra_hrs) + float(ra_min)/60 + float(ra_sec)/3600) * 15
    
    # Dec: degrees
    dec_decimal = abs(float(dec_deg)) + float(dec_min)/60 + float(dec_sec)/3600
    if dec_sign == '-':
        dec_decimal = -dec_decimal
    
    return ra_decimal, dec_decimal

def ra_dec_to_cartesian(ra_deg, dec_deg, radius=1.0):
    """Convert RA/Dec (in degrees) to Cartesian coordinates on a sphere."""
    ra_rad = np.radians(ra_deg)
    dec_rad = np.radians(dec_deg)
    
    x = radius * np.cos(dec_rad) * np.cos(ra_rad)
    y = radius * np.cos(dec_rad) * np.sin(ra_rad)
    z = radius * np.sin(dec_rad)
    
    return x, y, z

def convert_constellation_boundaries(input_file, output_file, radius=1.0):
    """
    Convert fixed-width constellation boundary file to clean CSV format.
    
    Args:
        input_file: Path to the original fixed-width format file
        output_file: Path for the output CSV file
        radius: Radius for the Cartesian coordinates (default 1.0 for unit sphere)
    """
    
    boundaries = []
    
    # Read the fixed-width format file
    with open(input_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            
            try:
                # Parse fixed-width format
                vertex1_key = line[0:3].strip()
                vertex2_key = line[4:7].strip()
                edge_type = line[8]  # M or P
                edge_dir = line[9]   # + or -
                
                # First vertex
                ra1_hrs = line[11:13].strip() or "0"
                ra1_min = line[14:16].strip() or "0"
                ra1_sec = line[17:19].strip() or "0"
                dec1_sign = line[20]
                dec1_deg = line[21:23].strip() or "0"
                dec1_min = line[24:26].strip() or "0"
                dec1_sec = line[27:29].strip() or "0"
                
                # Second vertex
                ra2_hrs = line[30:32].strip() or "0"
                ra2_min = line[33:35].strip() or "0"
                ra2_sec = line[36:38].strip() or "0"
                dec2_sign = line[39]
                dec2_deg = line[40:42].strip() or "0"
                dec2_min = line[43:45].strip() or "0"
                dec2_sec = line[46:48].strip() or "0"
                
                # Constellation
                constellation = line[49:57].strip()
                
                # Parse coordinates
                ra1, dec1 = parse_ra_dec(ra1_hrs, ra1_min, ra1_sec, dec1_sign, dec1_deg, dec1_min, dec1_sec)
                ra2, dec2 = parse_ra_dec(ra2_hrs, ra2_min, ra2_sec, dec2_sign, dec2_deg, dec2_min, dec2_sec)
                ra1 +=1.418
                ra2 +=1.418
                # Convert to Cartesian
                x1, y1, z1 = ra_dec_to_cartesian(ra1, dec1, radius)
                x2, y2, z2 = ra_dec_to_cartesian(ra2, dec2, radius)
                
                boundaries.append({
                    'vertex1_key': vertex1_key,
                    'vertex2_key': vertex2_key,
                    'edge_type': edge_type,
                    'edge_dir': edge_dir,
                    'ra1_deg': ra1,
                    'dec1_deg': dec1,
                    'ra2_deg': ra2,
                    'dec2_deg': dec2,
                    'x1': x1,
                    'y1': y1,
                    'z1': z1,
                    'x2': x2,
                    'y2': y2,
                    'z2': z2,
                    'constellation': constellation
                })
            except Exception as e:
                print(f"Error parsing line: {line}")
                print(f"Error: {e}")
                continue
    
    # Write to CSV
    with open(output_file, 'w', newline='') as f:
        fieldnames = ['vertex1_key', 'vertex2_key', 'edge_type', 'edge_dir',
                     'ra1_deg', 'dec1_deg', 'ra2_deg', 'dec2_deg',
                     'x1', 'y1', 'z1', 'x2', 'y2', 'z2', 'constellation']
        
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(boundaries)
    
    print(f"Converted {len(boundaries)} boundaries")
    print(f"Output saved to: {output_file}")
    return boundaries


# Example usage

# convert_constellation_boundaries(
#     "borders.txt",  # Your input file
#     "constellation_boundaries_clean.csv",     # Output file
#     radius=1.0  # Use 1.0 for unit sphere, scale later in Manim
# )

print("\n✓ Conversion complete!")
print("Now you can use the clean CSV with the simplified loader below.")

import pandas as pd

# Caminho do arquivo CSV
csv_path = "estrelas.csv"

# Lê o arquivo
df = pd.read_csv("data_set/constellation_boundaries_correct.csv")

# Ordena pelo valor de RA (ascensão reta)
df_sorted = df.sort_values(by="ra1_deg")

# Salva o resultado num novo arquivo
df_sorted.to_csv("estrelas_sorted.csv", index=False)

print("Arquivo ordenado salvo como 'estrelas_sorted.csv'.")