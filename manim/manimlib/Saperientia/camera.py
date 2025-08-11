import numpy as np
from manimlib.constants import *
from manimlib.utils.rate_functions import smooth


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
