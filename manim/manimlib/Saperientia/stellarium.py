from manimlib import *
import pandas as pd
import numpy as np
import math
from manimlib.Saperientia.Variables import * 
import csv


def rgb_to_hex(rgb):
    return '#' + ''.join(f'{int(round(c * 255)):02X}' for c in rgb)

def inverse_interpolate(x0, x1, x):
    return (x - x0) / (x1 - x0)


def interpolate_color(c1, c2, t):
    return [
        (1 - t) * a + t * b
        for a, b in zip(c1, c2)
    ]


def hex_to_rgb(hex_color):
    hex_color = hex_color.lstrip('#')
    return [int(hex_color[i:i+2], 16) / 255 for i in (0, 2, 4)]


def color_index_to_hex(bv_index):
    colors_hex = [
        "#E6F0FF",
        "#B3D1FF",
        "#80B2FF",
        "#FFFFFF",
        "#FFFF99",
        "#FFFF00",
        "#FFCC00",
        "#FF9900",
        "#FF6600",
        "#FF3300",
        "#FF0000",
    ]
    
    colors = [hex_to_rgb(c) for c in colors_hex]

    alpha = inverse_interpolate(-0.2, 2.9, bv_index)
    alpha = min(max(alpha, 0), 1)

    n_segments = len(colors) - 1
    scaled_alpha = alpha * n_segments
    idx = int(scaled_alpha)
    idx = min(idx, n_segments - 1)
    local_alpha = scaled_alpha - idx

    c1 = colors[idx]
    c2 = colors[idx + 1]

    rgb = interpolate_color(c1, c2, local_alpha)
    return rgb_to_hex(rgb)


def extract_star_data(file_path='data_set/sorted_hygdata_v41.csv',lat=0,lon=0):
    stars = []

    try:
        df = pd.read_csv(file_path, dtype={5: str, 27: str})
    except FileNotFoundError:
        print(f"File not found: {file_path}")
        return []

    for _, row in df.head(STAR_NUMBER).iterrows():
        try:
            dec = float(row['dec'])
            ra = float(row['ra']) * 15  # RA convertido de horas para graus
            mag = float(row['mag'])
            distance = float(row['dist'])
            ci = row['ci']
            lum = row['lum']
            
        except Exception as e:
            print(f"Error processing star data: {e}")
            continue

        if distance <= 0 :#or mag > 6
            continue

        # theta = np.pi / 2 - np.radians(dec)
        # phi = np.radians(ra) - np.pi / 2

        if not np.isnan(ci):
            color = color_index_to_hex(ci)
        else:
            color = "#ffffff"

        x = 3.08e16 * row['x'] / MANIM_REAL_SCALE
        y = 3.08e16 * row['y'] / MANIM_REAL_SCALE
        z = 3.08e16 * row['z'] / MANIM_REAL_SCALE
        
        

        # if row['proper'] == 'Alnilam':
        #     print(x, y, z,theta,phi)
        # x = 10 * math.sin(theta) * math.cos(phi)
        # y = 10 * math.sin(theta) * math.sin(phi)
        # z = 10 * math.cos(theta)
    
        # color = color_index_to_rgb(ci) if not np.isnan(ci) else [1, 1, 1]
        
        distance = 3.08e16 * distance / MANIM_REAL_SCALE 
        coeficiente = -0.4
        size = distance  * np.exp(coeficiente * mag) # Ajuste de escala para visualização VALOR
        stars.append({
            "x": x,
            "y": y,
            "z": z,
            "size": size,
            "color": color,
            "name": row["proper"],
            "hip": row["hip"] if "hip" in row else None,
            "mag": mag,
            'dec': dec,
            'ra': ra,
            "distance": distance,
            "ci": ci,
            "lum": lum,
            "vx": row["vx"],
            "vy": row["vy"],
            "vz": row["vz"],
        
        })

    return stars

def Asterisms(
            sphere_radius = 10000,
            opacity = 0.3,
            time_years = 0,
            constellations = None,
            stroke_width = None,
            lat=90,
            lon=0,
            ):
    
    if stroke_width is None:
        stroke_width = LINE_SIZE*sphere_radius
    
    stars_data = extract_star_data()
    star_lookup = {}

    # Update star positions with proper motion
    for star_dic in stars_data:
        x, y, z = star_dic["x"], star_dic["y"], star_dic["z"]
        vx, vy, vz = star_dic["vx"], star_dic["vy"], star_dic["vz"]
        x += time_years * vx * 3.08e16 / MANIM_REAL_SCALE
        y += time_years * vy * 3.08e16 / MANIM_REAL_SCALE
        z += time_years * vz * 3.08e16 / MANIM_REAL_SCALE
        hip = star_dic["hip"]
        star_dic["x"] = x
        star_dic["y"] = y
        star_dic["z"] = z
        star_lookup[hip] = star_dic
    
    if constellations:
        constellations = [c.upper() for c in constellations]
    
    lines = Group()
    
    # Process each constellation's asterism data
    for asterism in ASTERISM_DATA:
        const_code = asterism[0]
        
        # Filter by constellation if specified
        if constellations and const_code.upper() not in constellations:
            continue
        
        # Extract HIP numbers (skip first two elements: code and line count)
        hip_numbers = asterism[2:]
        
        # Create lines by connecting consecutive pairs
        for i in range(0, len(hip_numbers) - 1, 2):
            hip1 = hip_numbers[i]
            hip2 = hip_numbers[i + 1]
            
            if hip1 not in star_lookup or hip2 not in star_lookup:
                continue

            s1 = star_lookup[hip1]
            s2 = star_lookup[hip2]

            p1 = np.array([s1["x"], s1["y"], s1["z"]])
            p2 = np.array([s2["x"], s2["y"], s2["z"]])
            
            # Normalize to sphere radius
            norm1 = np.linalg.norm(p1)
            norm2 = np.linalg.norm(p2)
            
            if norm1 > 0 and norm2 > 0:
                p1 = p1 * sphere_radius / norm1
                p2 = p2 * sphere_radius / norm2

                line = Line(
                    p1, p2,
                    color=WHITE,
                    opacity=opacity,
                ).set_stroke(width=stroke_width,opacity=opacity).set_z_index(-5).deactivate_depth_test()
                lines.add(line)
    lines.rotate(lon*DEG,axis=OUT).rotate((lat-90)*DEG,axis=RIGHT,about_point=ORIGIN)
    return lines

#Testar com Arc
def Stars(  sphere_mode = False, 
            sphere_radius = 1,
            size_factor = 1,
            star_opacity = 1,
            ayre = False,
            glow = True,
            hr_mode = False,
            hr_mode_radius = False,
            time_years = 0,
            lat=90,
            lon=0,
            ):
    
    stars_data = extract_star_data()
    stars = Group()

    for star_dic in stars_data:
        x, y, z = star_dic["x"], star_dic["y"], star_dic["z"]
        size = star_dic["size"]
        color = star_dic["color"]
        # hip = star_dic["hip"]
        distance = star_dic["distance"]
        lum = star_dic["lum"]
        ci=star_dic["ci"]
        
        if sphere_mode:
            pos_factor=sphere_radius/distance
            size_pos_factor = size_factor/distance*sphere_radius*0.05
            opacity=star_opacity
        elif hr_mode:
            pos_factor=sphere_radius
            size_pos_factor = size_factor/distance*sphere_radius*0.05
            opacity=star_opacity
            z=0
            y=np.log10(lum)/2-1.5
            x=(ci-0.656)*2
        elif hr_mode_radius:
            pos_factor=sphere_radius
            size_pos_factor = size_factor/distance*sphere_radius*0.05
            opacity=star_opacity
            T = 4600*((1/(0.92*ci+1.7))+(1/(0.92*ci+0.62)))
            z=np.sqrt(lum*3.8e26/(4*PI*5.67e-8*T**4))/695_500_000/15
            y=np.log10(lum)/2-1.5
            x=(ci-0.656)*2
        elif ayre:
            pos_factor=sphere_radius/distance
            size_pos_factor = size_factor/distance*sphere_radius*0.05
            if z >= 0:
                xy = np.sqrt(x**2 + y**2)
                theta=np.arctan2(xy,z)
                x=x/xy*theta*distance
                y=y/xy*theta*distance
                z=sphere_radius
                opacity=star_opacity
            else:
                opacity=0
        else:
            pos_factor=1
            x+=time_years * star_dic["vx"]* 3.08e16 / MANIM_REAL_SCALE
            y+=time_years * star_dic["vy"]* 3.08e16 / MANIM_REAL_SCALE
            z+=time_years * star_dic["vz"]* 3.08e16 / MANIM_REAL_SCALE
            size_pos_factor = size_factor/70
            opacity=star_opacity
            
        if glow: 
            star = GlowDots(
                points=[np.array([x, y, z])*pos_factor],
                color=color,
                radius=size*size_pos_factor,
                opacity=opacity
            )
        else:
            star = DotCloud(
                points=[np.array([x, y, z])*pos_factor],
                color=color,
                radius=size*size_pos_factor/3,
                opacity=opacity
            )
        stars.add(star)
    stars.rotate(lon*DEG,axis=OUT).rotate((lat-90)*DEG,axis=RIGHT,about_point=ORIGIN)
    return stars


def load_constellation_boundaries_clean(csv_file="data_set/constellations.csv"):
    """Load pre-processed constellation boundaries from correct CSV."""
    boundaries = []
    
    with open(csv_file, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            
            boundaries.append({
                'vertex1_key': row['vertex1_key'],
                'vertex2_key': row['vertex2_key'],
                'edge_type': row['edge_type'],
                'edge_dir': row['edge_dir'],
                'ra1': float(row['ra1_deg']),
                'dec1': float(row['dec1_deg']),
                'ra2': float(row['ra2_deg']),
                'dec2': float(row['dec2_deg']),
                'x1': float(row['x1']),
                'y1': float(row['y1']),
                'z1': float(row['z1']),
                'x2': float(row['x2']),
                'y2': float(row['y2']),
                'z2': float(row['z2']),
                'constellation': row['constellation']
            })
    
    return boundaries

def create_boundary_edge(boundary, radius, num_points_dec=2,num_points_ra=8):
    """Create points along a boundary edge (meridian or parallel)."""
    edge_type = boundary['edge_type']
    
    # Get RA/Dec for interpolation
    ra1, dec1 = boundary['ra1'], boundary['dec1']
    ra2, dec2 = boundary['ra2'], boundary['dec2']
    
    points = []
    
    if edge_type == 'M':  # Meridian (constant RA)
        ra = ra1
        decs = np.linspace(dec1, dec2, num_points_dec)
        for dec in decs:
            ra_rad = np.radians(ra)
            dec_rad = np.radians(dec)
            x = radius * np.cos(dec_rad) * np.cos(ra_rad)
            y = radius * np.cos(dec_rad) * np.sin(ra_rad)
            z = radius * np.sin(dec_rad)
            points.append(np.array([x, y, z]))
    
    elif edge_type == 'P':  # Parallel (constant Dec)
        dec = dec1
        # Handle RA wrapping around 360/0 boundary
        ra_start = ra1
        ra_end = ra2
        
        # Check if we need to wrap around
        if abs(ra_end - ra_start) > 180:
            if ra_start < ra_end:
                ra_start += 360
            else:
                ra_end += 360
        
        ras = np.linspace(ra_start, ra_end, num_points_ra)
        dec_rad = np.radians(dec)
        
        for ra in ras:
            ra_wrapped = ra % 360
            ra_rad = np.radians(ra_wrapped)
            x = radius * np.cos(dec_rad) * np.cos(ra_rad)
            y = radius * np.cos(dec_rad) * np.sin(ra_rad)
            z = radius * np.sin(dec_rad)
            points.append(np.array([x, y, z]))
    
    return points

def Constellations(csv_file="data_set/constellations.csv", sphere_radius=10000, constellations=None, 
                   color=BLUE, stroke_width=None, verbose=False,lat=90,lon=0):
    """
    Fast function to create constellation boundary lines from pre-processed CSV.
    
    Args:
        csv_file: Path to the correct CSV file (pre-processed)
        radius: Radius of the sphere
        constellations: Single constellation abbreviation (e.g., "AND"), 
                       list of constellations (e.g., ["AND", "ORI", "UMA"]),
                       or None for all
        color: Line color (default BLUE)
        stroke_width: Line width (default 2)
        verbose: Print debug information
    
    Returns:
        Group of Line3D objects
    """
    all_boundaries = load_constellation_boundaries_clean(csv_file)
    if stroke_width is None:
        stroke_width = LINE_SIZE*sphere_radius
    
    # Filter by constellation(s) if specified
    if constellations:
        # Convert single string to list for uniform handling
        if isinstance(constellations, str):
            constellations = [constellations]
        
        # Convert to uppercase for case-insensitive matching
        constellations = [c.upper() for c in constellations]
        
        # Filter boundaries
        boundaries = [b for b in all_boundaries 
                     if any(const in b['constellation'].upper() for const in constellations)]
        
        if verbose:
            print(f"Found {len(boundaries)} edges for {constellations}")
    else:
        boundaries = all_boundaries
        if verbose:
            print(f"Loaded {len(boundaries)} total edges")
    
    if verbose:
        constellations_found = set(b['constellation'] for b in boundaries)
        print(f"Constellations in filtered set: {sorted(constellations_found)}")
    
    lines = Group()
    
    for boundary in boundaries:
        points = create_boundary_edge(boundary, sphere_radius)
        
        if len(points) < 2:
            continue
        
        for i in range(len(points) - 1):
            line = Line(
                points[i],
                points[i + 1],
                color=color,
            ).set_stroke(width=stroke_width)
            lines.add(line)
    
    if verbose:
        print(f"Created {len(lines)} line segments")
    lines.rotate(lon*DEG,axis=OUT).rotate((lat-90)*DEG,axis=RIGHT,about_point=ORIGIN)
    return lines
