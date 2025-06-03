from manimlib import *
from manimlib.Saperientia.astro_structures import * 
from manimlib.Saperientia.stellarium import * 
import numpy as np
from scipy.spatial.transform import Rotation as R

class Teste(ThreeDScene):
    def construct(self):
        self.frame.reorient(0, 0, 0, ORIGIN+[0,0,0.001], 0.0005)
        # self.camera.frame.set_width(10000)
        stars_data = extract_star_data()

        stars = Group()

        for x, y, z, size, color in stars_data:
            star = DotCloud(
                points=[np.array([x, y, z])],
                color=color,
                radius=size,
                opacity=0.75
            )
            stars.add(star)

        self.add(stars)
    
        sun = Sun()
        sun.add_updater(lambda m:self.add(sun))
        earth = Earth().shift(RIGHT *20)
        light = self.camera.light_source
        moon = Moon()  # Lua ainda menor

        self.add(sun, earth,moon)
        
        

        # Trackers para os ângulos da Terra e da Lua
        angle_earth = ValueTracker(0)
        angle_moon = ValueTracker(0)

        # Atualiza a posição da Terra orbitando o Sol
        earth_orbit = always_redraw(lambda: 
            earth.move_to(sun.get_center() + radius_orbit_earth * np.array([
                np.cos(angle_earth.get_value()),
                np.sin(angle_earth.get_value()),
                0
            ]))
        )
        
        # Atualiza a posição da Lua orbitando a Terra
        moon_orbit = always_redraw(lambda:
            moon.move_to(earth.get_center() + radius_orbit_moon * np.array([
                np.cos(angle_moon.get_value()),
                np.sin(angle_moon.get_value()),
                0
            ]))
        )
        # self.frame.reorient(330,70,0,[53.38888931,  0.        ,  0.        ],15)
        self.frame.reorient(336,70,0,[53.38888931,  0.        ,  0.        ],15)
        esfera_celeste = EsferaCeleste().move_to([53.38888931,  0.        ,  0.        ])
        equador = Equador(90).move_to([53.38888931,  0.        ,  0.        ])
        eixopolar = EixoPolar(90).move_to([53.38888931,  0.        ,  0.        ])
        esfera_celeste.add_updater(lambda m:self.add(esfera_celeste))
        self.play(ShowCreation(esfera_celeste),run_time=2)
        self.play(ShowCreation(equador),ShowCreation(eixopolar),run_time=4)
        self.play(Rotate(earth, 2*TAU,axis=Z_AXIS),run_time=4)
        p1 = P(-30,200,center=[53.38888931,  0.        ,  0.        ],tamanho=2*DEFAULT_DOT_RADIUS,cor=RED)
        p2 = P(30,190,center=[53.38888931,  0.        ,  0.        ],tamanho=2*DEFAULT_DOT_RADIUS,cor=GREEN)
        p3 = P(10,240,center=[53.38888931,  0.        ,  0.        ],tamanho=2*DEFAULT_DOT_RADIUS,cor=YELLOW)
        self.play(FadeIn(p1),FadeIn(p2),FadeIn(p3))
        arc1 = GrandeArco(p1,p2,cor=PINK, espessura=4,center=[53.38888931,  0.        ,  0.        ])
        arc2 = GrandeArco(p1,p3,cor=PINK, espessura=4,center=[53.38888931,  0.        ,  0.        ])
        arc3 = GrandeArco(p2,p3,cor=PINK, espessura=4,center=[53.38888931,  0.        ,  0.        ]) 
        
        self.play(ShowCreation(arc1),ShowCreation(arc2),ShowCreation(arc3))
        self.play(self.frame.animate.reorient(330,110,0,[53.38888931,  0.        ,  0.        ],15),run_time=2)
        self.play(self.frame.animate.reorient(336,110,0,[53.38888931,  0.        ,  0.        ],15),run_time=2)
        self.wait(2)
        self.play(self.frame.animate.reorient(330,70,0,[53.38888931,  0.        ,  0.        ],15),run_time=2)
        ang1 = AnguloEsferico(p1,p3,p2,espessura=3,raio_circulo=0.7,center=[53.38888931,  0.        ,  0.        ])
        ang2 = AnguloEsferico(p3,p2,p1,espessura=3,raio_circulo=0.7,center=[53.38888931,  0.        ,  0.        ])
        ang3 = AnguloEsferico(p2,p1,p3,espessura=3,raio_circulo=0.7,center=[53.38888931,  0.        ,  0.        ])
        self.play(ShowCreation(ang1),ShowCreation(ang2),ShowCreation(ang3))
        self.play(
            angle_moon.animate.increment_value(4 * PI),
            self.frame.animate.reorient(60,80,0,[53.38888931,  0.        ,  0.        ],30),
            run_time=5,
        )
        self.play(FadeOut(arc1),FadeOut(arc2),FadeOut(arc3),FadeOut(ang1),FadeOut(ang2),FadeOut(ang3),FadeOut(equador),FadeOut(p1),FadeOut(p2),FadeOut(p3),FadeOut(eixopolar),FadeOut(esfera_celeste))
        self.play(
            self.frame.animate.reorient(90,80,0,[0,  0.        ,  0.        ],90),
            angle_earth.animate.increment_value(2 * PI),
            angle_moon.animate.increment_value(8 * 2 * PI),
            Rotate(earth,10*TAU,axis=Z_AXIS),
            run_time=20,
            rate_func=linear
        )
        
        self.play(
            self.frame.animate.increment_theta(-1.5*PI),
            angle_earth.animate.increment_value(2 * PI),
            angle_moon.animate.increment_value(8 * 2 * PI),
            Rotate(earth,10*TAU,axis=Z_AXIS),
            run_time=20,
            rate_func=linear
        )
        
        self.play(self.frame.animate.reorient(10,80,0,[0,  0.        ,  0.        ],300000),run_time=10)
        self.play(self.frame.animate.reorient(10,80,0,[300000,  300000.        ,  100.        ],300000),run_time=5)
        
        
        
        
        # p1 = [ 35.85023998, -29.6991819 ,   2.47037005]
        
        # self.play(
        #     angle_earth.animate.increment_value(2 * PI),
        #     angle_moon.animate.increment_value(8 * 2 * PI),
        #     run_time=20,
        #     rate_func=linear
        # )
        self.embed()
        
class Olho_direito(ThreeDScene):
    def construct(self):
        self.frame.reorient(0, 0, 0, ORIGIN+[0,0,0.001], 0.0005)
        # self.camera.frame.set_width(10000)
        stars_data = extract_star_data()

        stars = Group()

        for x, y, z, size, color in stars_data:
            star = DotCloud(
                points=[np.array([x, y, z])],
                color=color,
                radius=size,
                opacity=0.75
            )
            stars.add(star)

        self.add(stars)
    
        sun = Sun()
        sun.add_updater(lambda m:self.add(sun))
        earth = Earth().shift(RIGHT *20)
        light = self.camera.light_source
        moon = Moon()  # Lua ainda menor

        self.add(sun, earth,moon)
        
        

        # Trackers para os ângulos da Terra e da Lua
        angle_earth = ValueTracker(0)
        angle_moon = ValueTracker(0)

        # Atualiza a posição da Terra orbitando o Sol
        earth_orbit = always_redraw(lambda: 
            earth.move_to(sun.get_center() + radius_orbit_earth * np.array([
                np.cos(angle_earth.get_value()),
                np.sin(angle_earth.get_value()),
                0
            ]))
        )
        
        # Atualiza a posição da Lua orbitando a Terra
        moon_orbit = always_redraw(lambda:
            moon.move_to(earth.get_center() + radius_orbit_moon * np.array([
                np.cos(angle_moon.get_value()),
                np.sin(angle_moon.get_value()),
                0
            ]))
        )
        # self.frame.reorient(330,70,0,[53.38888931,  0.        ,  0.        ],15)
        self.frame.reorient(336,70,0,[53.38888931,  0.        ,  0.        ],15)
        esfera_celeste = EsferaCeleste().move_to([53.38888931,  0.        ,  0.        ])
        equador = Equador(90).move_to([53.38888931,  0.        ,  0.        ])
        eixopolar = EixoPolar(90).move_to([53.38888931,  0.        ,  0.        ])
        esfera_celeste.add_updater(lambda m:self.add(esfera_celeste))
        self.play(ShowCreation(esfera_celeste),run_time=2)
        self.play(ShowCreation(equador),ShowCreation(eixopolar),run_time=4)
        self.play(Rotate(earth, 2*TAU,axis=Z_AXIS),run_time=4)
        p1 = P(-30,200,center=[53.38888931,  0.        ,  0.        ],tamanho=2*DEFAULT_DOT_RADIUS,cor=RED)
        p2 = P(30,190,center=[53.38888931,  0.        ,  0.        ],tamanho=2*DEFAULT_DOT_RADIUS,cor=GREEN)
        p3 = P(10,240,center=[53.38888931,  0.        ,  0.        ],tamanho=2*DEFAULT_DOT_RADIUS,cor=YELLOW)
        self.play(FadeIn(p1),FadeIn(p2),FadeIn(p3))
        arc1 = GrandeArco(p1,p2,cor=PINK, espessura=4,center=[53.38888931,  0.        ,  0.        ])
        arc2 = GrandeArco(p1,p3,cor=PINK, espessura=4,center=[53.38888931,  0.        ,  0.        ])
        arc3 = GrandeArco(p2,p3,cor=PINK, espessura=4,center=[53.38888931,  0.        ,  0.        ]) 
        
        self.play(ShowCreation(arc1),ShowCreation(arc2),ShowCreation(arc3))
        self.play(self.frame.animate.reorient(336,110,0,[53.38888931,  0.        ,  0.        ],15),run_time=2)
        self.wait(2)
        self.play(self.frame.animate.reorient(336,70,0,[53.38888931,  0.        ,  0.        ],15),run_time=2)
        ang1 = AnguloEsferico(p1,p3,p2,espessura=3,raio_circulo=0.7,center=[53.38888931,  0.        ,  0.        ])
        ang2 = AnguloEsferico(p3,p2,p1,espessura=3,raio_circulo=0.7,center=[53.38888931,  0.        ,  0.        ])
        ang3 = AnguloEsferico(p2,p1,p3,espessura=3,raio_circulo=0.7,center=[53.38888931,  0.        ,  0.        ])
        self.play(ShowCreation(ang1),ShowCreation(ang2),ShowCreation(ang3))
        self.play(
            angle_moon.animate.increment_value(4 * PI),
            self.frame.animate.reorient(63,80,0,[53.38888931,  0.        ,  0.        ],30),
            run_time=5,
        )
        self.play(FadeOut(arc1),FadeOut(arc2),FadeOut(arc3),FadeOut(ang1),FadeOut(ang2),FadeOut(ang3),FadeOut(equador),FadeOut(p1),FadeOut(p2),FadeOut(p3),FadeOut(eixopolar),FadeOut(esfera_celeste))
        self.play(
            self.frame.animate.reorient(91,80,0,[0,  0.        ,  0.        ],90),
            angle_earth.animate.increment_value(2 * PI),
            angle_moon.animate.increment_value(8 * 2 * PI),
            Rotate(earth,10*TAU,axis=Z_AXIS),
            run_time=20,
            rate_func=linear
        )
        
        self.play(
            self.frame.animate.increment_theta(-1.5*PI),
            angle_earth.animate.increment_value(2 * PI),
            angle_moon.animate.increment_value(8 * 2 * PI),
            Rotate(earth,10*TAU,axis=Z_AXIS),
            run_time=20,
            rate_func=linear
        )
        
        self.play(self.frame.animate.reorient(10.0003,80,0,[0,  0.        ,  0.        ],300000),run_time=10)
        self.play(self.frame.animate.reorient(10.0003,80,0,[300000,  300000.        ,  100.        ],300000),run_time=5)
        
def convert_camera_angles(phi_p, theta_p, gamma_p, alpha, degrees=True):
    if degrees:
        phi_p = np.radians(phi_p)
        theta_p = np.radians(theta_p)
        gamma_p = np.radians(gamma_p)
        alpha = np.radians(alpha)
    
    # Rotação da câmera na base inclinada (sistema ZXZ)
    R_p = R.from_euler('ZXZ', [theta_p, phi_p, gamma_p])
    
    # Rotação da base inclinada em relação à base original
    R_tilt = R.from_euler('X', alpha)

    # Compor as rotações (usando * ao invés de @)
    R_total = R_tilt * R_p

    # Obter ângulos na base original
    theta, phi, gamma = R_total.as_euler('ZXZ')

    if degrees:
        return np.degrees(phi), np.degrees(theta), np.degrees(gamma)
    else:
        return phi, theta, gamma

class Terra(ThreeDScene):
    def construct(self):
        self.frame.reorient(0,130,0,ORIGIN,15)
        nuvem = Clouds().set_opacity(0.2)
        terra = Earth()
        esfera = EsferaCeleste().set_shading(1,1,1)
        self.add(terra)
        # self.play(self.frame.animate.reorient(180,90,0,[0,  0       ,  0.        ],5),run_time=5)
        # self.play(self.frame.animate.reorient(91,90,-90,[0,  3.        ,  0.        ],0.01),run_time=5)
        p1 = P(-23,180,raio=3).get_center()
        stars_data = extract_star_data()

        stars = Group()

        for x, y, z, size, color in stars_data:
            star = DotCloud(
                points=[np.array([x, y, z])],
                color=color,
                radius=size,
                opacity=0.75
            )
            stars.add(star)
        self.add(stars)
        
        terra.rotate_about_origin(150*DEGREES,Z_AXIS)
        self.frame.set_field_of_view(1.3)
        self.wait(3)
        phi,theta,gamma=convert_camera_angles(80,180,0,90+23)
        self.play(self.frame.animate.reorient(theta,phi,gamma,p1*1.07,0.5),run_time=3)
        phi,theta,gamma=convert_camera_angles(110,180,0,90+23)
        self.play(self.frame.animate.reorient(theta,phi,gamma,p1*1.001,0.001),run_time=3)
        self.play(stars.animate.rotate(60*DEGREES,axis=Z_AXIS, about_point=ORIGIN),rate_func=linear)
        phi,theta,gamma=convert_camera_angles(110,0,0,90+23)
        self.play(self.frame.animate.reorient(theta,phi,gamma,p1*1.001,0.001),run_time=3)
        phi,theta,gamma=convert_camera_angles(150,0,0,90+23)
        self.play(self.frame.animate.reorient(theta,phi,gamma,p1*1.001,0.001),run_time=3)
        phi,theta,gamma=convert_camera_angles(110,0,0,90+23)
        self.play(self.frame.animate.reorient(theta,phi,gamma,p1*1.001,0.001),run_time=3)
        phi,theta,gamma=convert_camera_angles(110,270,0,90+23)
        self.play(self.frame.animate.reorient(theta,phi,gamma,p1*1.001,0.001),run_time=3)
        
        phi,theta,gamma=convert_camera_angles(80,180,0,90+23)
        self.play(self.frame.animate.reorient(theta,phi,gamma,p1*1.07,0.5),run_time=3)
        self.play(self.frame.animate.reorient(0,130,0,ORIGIN,10),run_time=3)
        
        
        camera = self.frame
        # self.play(self.frame.animate.reorient(180,90,0,p1,10),run_time=5)
        # self.play(self.frame.animate.reorient(120,90,-113,p1*1.0001,0.000001),run_time=5,rate_func=linear)
        # self.play(self.frame.animate.reorient(120,270,-67,p1*1.0001,0.000001),run_time=5,rate_func=linear)
        self.embed()
