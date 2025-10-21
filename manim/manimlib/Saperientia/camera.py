import numpy as np
from manimlib.utils.rate_functions import smooth
from manimlib import *
from manimlib.Saperientia.Variables import *
from manimlib.Saperientia.stellarium import *

def reorientar_camera(frame, 
                          point,
                          theta_deg=0, phi_deg=0, gamma_deg=0,
                          epsilon=1e-7,
                          fovy=45):
    """Inverte theta e phi no esquema zxz e cola a câmera no ponto usando epsilon minúsculo."""
    axes = getattr(frame, "euler_axes", "zxz")
    if axes != "zxz":
        raise ValueError(f"Sequência de Euler {axes} não implementada aqui")

    # Inversão correta
    th = theta_deg +180
    ph = 180.0 - phi_deg
    ga = gamma_deg

    # Reorienta com centro no ponto
    frame.reorient(theta_degrees=th, phi_degrees=ph, gamma_degrees=ga,center=np.array(point, dtype=float),height=epsilon)

    frame.uniforms["fovy"]=fovy*DEG

from manimlib.animation.animation import Animation

def _shortest_lerp_deg(a, b, t):
    d = (b - a + 180.0) % 360.0 - 180.0
    return a + t * d

class AnimarCamera(Animation):
    def __init__(self, frame, 
                 point,
                 theta_deg=None, phi_deg=None, gamma_deg=None,
                 epsilon=1e-6, fovy=None,
                 run_time=2, rate_func=smooth):
        self.frame = frame
        self.point_target = np.array(point, dtype=float)

        # estado atual (em graus)
        th0, ph0, ga0 = frame.get_euler_angles() * 180.0 / PI
        self.center0  = frame.get_center().copy()
        self.height0  = float(frame.get_height())
        self.fovy0    = float(frame.get_field_of_view() / DEG)

        # converte para espaço pré flip
        self.theta_pre0 = th0 - 180.0
        self.phi_pre0   = 180.0 - ph0
        self.gamma_pre0 = ga0

        # define alvos: se None, mantém o valor atual
        self.theta_t = float(theta_deg) if theta_deg is not None else self.theta_pre0
        self.phi_t   = float(phi_deg)   if phi_deg   is not None else self.phi_pre0
        self.gamma_t = float(gamma_deg) if gamma_deg is not None else self.gamma_pre0

        self.eps_t   = float(epsilon)
        self.fovy_t  = float(fovy) if fovy is not None else self.fovy0

        super().__init__(frame, run_time=run_time, rate_func=rate_func)

    def interpolate_mobject(self, alpha):
        th_pre = _shortest_lerp_deg(self.theta_pre0, self.theta_t, alpha)
        ph_pre = _shortest_lerp_deg(self.phi_pre0,   self.phi_t,   alpha)
        ga_pre = _shortest_lerp_deg(self.gamma_pre0, self.gamma_t, alpha)

        center = (1 - alpha) * self.center0 + alpha * self.point_target
        h      = (1 - alpha) * self.height0 + alpha * self.eps_t
        fv     = (1 - alpha) * self.fovy0   + alpha * self.fovy_t

        reorientar_camera(
            self.frame,
            point=center,
            theta_deg=th_pre, phi_deg=ph_pre, gamma_deg=ga_pre,
            epsilon=max(h, 1e-9),
            fovy=fv
        )


def on_mouse_scroll(
    self,
    point,
    offset,
    x_pixel_offset: float,
    y_pixel_offset: float
) -> None:
    event_data = {"point": point, "offset": offset}
    propagate_event = EVENT_DISPATCHER.dispatch(EventType.MouseScrollEvent, **event_data)
    if propagate_event is not None and propagate_event is False:
        return

    rel_offset = y_pixel_offset / self.camera.get_pixel_height()
    self.frame.set_field_of_view(self.frame.get_field_of_view() - self.scroll_sensitivity * rel_offset)

Scene.on_mouse_scroll = on_mouse_scroll


def look_at_coord(self, ra_deg, dec_deg,gamma=0,center=ORIGIN,height=1):
    """
    Points the camera at a specific celestial coordinate (RA, Dec).
    
    Args:
        ra_deg: Right Ascension in degrees (0-360)
        dec_deg: Declination in degrees (-90 to +90)
    """
    self.frame.reorient(ra_deg - 90, dec_deg + 90, gamma, center, height)

Scene.look_at_coord = look_at_coord

def look_at_star(self, hip_id, gamma=0, center=ORIGIN, height=1):
    """
    Points the camera at a specific star given its HIP (Hipparcos) catalog number.
    
    Args:
        hip_id: Hipparcos catalog number (integer)
        gamma: Roll angle in degrees (default 0)
        center: Center point (default ORIGIN)
        height: Camera height (default 1)
    """
    stars_data = extract_star_data()
    
    # Find the star with the given HIP ID
    star = None
    for star_dic in stars_data:
        if star_dic["hip"] == hip_id:
            star = star_dic
            break
    
    if star is None:
        print(f"Warning: Star with HIP {hip_id} not found")
        return
    
    # Get RA and Dec from the star data
    ra_deg = star["ra"]
    dec_deg = star["dec"]
    
    # Use look_at_coord to point at the star
    self.look_at_coord(ra_deg, dec_deg, gamma, center, height)

Scene.look_at_star = look_at_star

def look_at_constellation(self, constellation_code, gamma=0, center=ORIGIN, height=1):
    """
    Points the camera at a constellation by calculating its center from asterism stars.
    
    Args:
        constellation_code: 3-letter constellation code (e.g., "Aql", "Ori", "UMa")
        gamma: Roll angle in degrees (default 0)
        center: Center point (default ORIGIN)
        height: Camera height (default 1)
    """
    # Constellation data with their asterism lines
    # Find the constellation (case-insensitive)
    constellation_code = constellation_code.upper()
    asterism = None
    for ast in ASTERISM_DATA:
        if ast[0].upper() == constellation_code:
            asterism = ast
            break
    
    if asterism is None:
        print(f"Warning: Constellation '{constellation_code}' not found")
        return
    
    # Get all unique HIP IDs from the asterism
    hip_ids = set(asterism[2:])
    
    # Load star data
    stars_data = extract_star_data()
    star_lookup = {star["hip"]: star for star in stars_data}
    
    # Calculate average RA and Dec
    ra_sum = 0
    dec_sum = 0
    count = 0
    
    for hip_id in hip_ids:
        if hip_id in star_lookup:
            star = star_lookup[hip_id]
            ra_sum += star["ra"]
            dec_sum += star["dec"]
            count += 1
    
    if count == 0:
        print(f"Warning: No stars found for constellation '{constellation_code}'")
        return
    
    # Calculate center
    ra_center = ra_sum / count
    dec_center = dec_sum / count
    
    # Use look_at_coord to point at the constellation center
    self.look_at_coord(ra_center, dec_center, gamma, center, height)

Scene.look_at_constellation = look_at_constellation