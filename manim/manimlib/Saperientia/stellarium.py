from manimlib import *
import pandas as pd
import numpy as np
import math

RADIUS = 10

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


def extract_star_data(file_path='sorted_hygdata_v41.csv'):
    stars = []

    try:
        df = pd.read_csv(file_path)
    except FileNotFoundError:
        print(f"File not found: {file_path}")
        return []

    for _, row in df.head(6000).iterrows():
        try:
            dec = float(row['dec'])
            ra = float(row['ra']) * 15  # RA convertido de horas para graus
            mag = float(row['mag'])
            distance = float(row['dist'])
            ci = row['ci']
        except Exception as e:
            print(f"Error processing star data: {e}")
            continue

        if distance <= 0 or mag > 6:
            continue

        theta = np.pi / 2 - np.radians(dec)
        phi = np.radians(ra) - np.pi / 2

        if not np.isnan(ci):
            color = color_index_to_hex(ci)
        else:
            color = "#ffffff"

        x = 10000 * row['x']
        y = 10000 * row['y']
        z = 10000 * row['z']

        if row['proper'] == 'Altair':
            print(x, y, z)
        # x = 10 * math.sin(theta) * math.cos(phi)
        # y = 10 * math.sin(theta) * math.sin(phi)
        # z = 10 * math.cos(theta)
    
        # color = color_index_to_rgb(ci) if not np.isnan(ci) else [1, 1, 1]

        # size = 10 * np.exp(-0.11 * mag)  # Ajuste de escala para visualização
        size = distance * 70 * np.exp(-0.4 * mag) # Ajuste de escala para visualização

        stars.append((x, y, z, size, color))
    return stars


class StarField3D(ThreeDScene):
    def construct(self):
        self.frame.reorient(0, 0, 0, ORIGIN+[0,0,0.001], 0.0005)
        # self.camera.frame.set_width(10000)
        stars_data = extract_star_data()
        print('terminou')

        stars = Group()

        for x, y, z, size, color in stars_data:
            star = DotCloud(
                points=[np.array([x, y, z])],
                color=color,
                radius=size,
            )
            stars.add(star)

        self.add(stars)

        # def on_mouse_scroll(
        #     point,
        #     offset,
        #     x_pixel_offset: float,
        #     y_pixel_offset: float
        # ) -> None:
        #     event_data = {"point": point, "offset": offset}
        #     propagate_event = EVENT_DISPATCHER.dispatch(EventType.MouseScrollEvent, **event_data)
        #     if propagate_event is not None and propagate_event is False:
        #         return
        #     print(self.frame.get_field_of_view())
        #     # rel_offset = y_pixel_offset / self.camera.get_pixel_height()
        #     self.frame.set_field_of_view(
        #         max(0.5, min(1.2, self.frame.get_field_of_view() - y_pixel_offset / 50))
        #     )

        # self.frame.set_field_of_view(0.5)
        # self.on_mouse_scroll = on_mouse_scroll
        # self.frame.reorient(0, 0, 0, ORIGIN+[0,0,0.1], 0.05)
        # https://github.com/3b1b/videos/blob/e5a041d2094ca8f11e0cabd20e9dd199b581e3f5/_2025/cosmic_distance/planets.py
        # frame.reorient(0, 90, 0, ORIGIN, 3.42)
        
        # Posicionamento da câmera
        # self.frame.set_euler_angles(
        #     theta=0 * DEGREES,
        #     phi=70 * DEGREES
        # )
        self.play(self.frame.animate.reorient(90, 180, 0, ORIGIN+[0,0,0.001], 0.0005), run_time=5)
        # self.frame.set_euler_angles(phi=90, theta=0)

        self.embed()
        # self.wait(10)
