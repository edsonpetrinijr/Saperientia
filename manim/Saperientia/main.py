from manimlib import *
from manimlib.Saperientia.astro_structures import * 
from manimlib.Saperientia.stellarium import * 
import numpy as np
from scipy.spatial.transform import Rotation as R

class Teste(InteractiveScene):
    def construct(self):
        self.frame.reorient(0, 0, 0, ORIGIN+[0,0,0.001], 0.0005)
        # self.camera.frame.set_width(10000)
        # stars_data = extract_star_data()

        # stars = Group()

        # for x, y, z, size, color in stars_data:
        #     star = DotCloud(
        #         points=[np.array([x, y, z])],
        #         color=color,
        #         radius=size,
        #         opacity=0.75
        #     )
        #     stars.add(star)

        # self.add(stars)
    
        sun = Sun()
        #sun.add_updater(lambda m:self.add(sun))
        earth = Earth().shift(RIGHT * 20)
        light = self.camera.light_source
        moon = Moon()  # Lua ainda menor
        esfera_celeste = EsferaCeleste(opacidade=1).move_to(earth.get_center())
        self.play(ShowCreation(esfera_celeste),run_time=2)
        esfera_celeste.add_updater(lambda m:self.add(esfera_celeste))
        p1 = P(-30,200,center=[53.38888931, 0., 0.],tamanho=2*DEFAULT_DOT_RADIUS,cor=RED)
        # text.add_updater(lambda _:self.add(text))
        self.add(earth,moon)
        
        text = Tex("Ola, Edson",font_size=100).next_to(p1,RIGHT*3).shift(10*DOWN)
        text.rotate(90*DEGREES,X_AXIS)
        text.set_anti_alias_width(0)
        text.set_z_index(1)
        # text.add_updater(lambda m:self.add())
        self.add(text)

        # Trackers para os ângulos da Terra e da Lua
        angle_earth = ValueTracker(0)
        angle_moon = ValueTracker(0)

        # Atualiza a posição da Terra orbitando o Sol
        # earth_orbit = always_redraw(lambda: 
        #     earth.move_to(sun.get_center() + radius_orbit_earth * np.array([
        #         np.cos(angle_earth.get_value()),
        #         np.sin(angle_earth.get_value()),
        #         0
        #     ]))
        # )
        
        # # Atualiza a posição da Lua orbitando a Terra
        # moon_orbit = always_redraw(lambda:
        #     moon.move_to(earth.get_center() + radius_orbit_moon * np.array([
        #         np.cos(angle_moon.get_value()),
        #         np.sin(angle_moon.get_value()),
        #         0
        #     ]))
        # )
        # self.frame.reorient(330,70,0,earth.get_center(),15)
        # self.frame.reorient(336,70,0,earth.get_center(),15)
        self.frame.reorient(336,70,0,earth.get_center(),15)
        equador = Equador(90).move_to(earth.get_center())
        eixopolar = EixoPolar(90).move_to(earth.get_center())
        self.play(ShowCreation(equador),ShowCreation(eixopolar),run_time=4)
        self.play(Rotate(earth, 2*TAU,axis=Z_AXIS),run_time=4)
        p2 = P(30,190,center=earth.get_center(),tamanho=2*DEFAULT_DOT_RADIUS,cor=GREEN)
        p3 = P(10,240,center=earth.get_center(),tamanho=2*DEFAULT_DOT_RADIUS,cor=YELLOW)
        self.play(FadeIn(p1),FadeIn(p2),FadeIn(p3))
        arc1 = GrandeArco(p1,p2,cor=PINK, espessura=4,center=earth.get_center())
        arc2 = GrandeArco(p1,p3,cor=PINK, espessura=4,center=earth.get_center())
        arc3 = GrandeArco(p2,p3,cor=PINK, espessura=4,center=earth.get_center())
        
        self.play(ShowCreation(arc1),ShowCreation(arc2),ShowCreation(arc3))

        # self.play(self.frame.animate.reorient(330,110,0,earth.get_center(),15),run_time=2)
        # self.play(self.frame.animate.reorient(336,110,0,earth.get_center(),15),run_time=2)
        # self.wait(2)
        # self.play(self.frame.animate.reorient(330,70,0,earth.get_center(),15),run_time=2)
        # ang1 = AnguloEsferico(p1,p3,p2,espessura=3,raio_circulo=0.7,center=earth.get_center())
        # ang2 = AnguloEsferico(p3,p2,p1,espessura=3,raio_circulo=0.7,center=earth.get_center())
        # ang3 = AnguloEsferico(p2,p1,p3,espessura=3,raio_circulo=0.7,center=earth.get_center())
        # self.play(ShowCreation(ang1),ShowCreation(ang2),ShowCreation(ang3))
        # self.play(
        #     angle_moon.animate.increment_value(4 * PI),
        #     self.frame.animate.reorient(60,80,0,earth.get_center(),30),
        #     run_time=5,
        # )
        # self.play(FadeOut(arc1),FadeOut(arc2),FadeOut(arc3),FadeOut(ang1),FadeOut(ang2),FadeOut(ang3),FadeOut(equador),FadeOut(p1),FadeOut(p2),FadeOut(p3),FadeOut(eixopolar),FadeOut(esfera_celeste))
        # self.play(
        #     self.frame.animate.reorient(90,80,0,[0,  0.        ,  0.        ],90),
        #     angle_earth.animate.increment_value(2 * PI),
        #     angle_moon.animate.increment_value(8 * 2 * PI),
        #     Rotate(earth,10*TAU,axis=Z_AXIS),
        #     run_time=20,
        #     rate_func=linear
        # )
        
        # self.play(
        #     self.frame.animate.increment_theta(-1.5*PI),
        #     angle_earth.animate.increment_value(2 * PI),
        #     angle_moon.animate.increment_value(8 * 2 * PI),
        #     Rotate(earth,10*TAU,axis=Z_AXIS),
        #     run_time=20,
        #     rate_func=linear
        # )
        
        # self.play(self.frame.animate.reorient(10,80,0,[0,  0.        ,  0.        ],300000),run_time=10)
        # self.play(self.frame.animate.reorient(10,80,0,[300000,  300000.        ,  100.        ],300000),run_time=5)

        # # p1 = [ 35.85023998, -29.6991819 ,   2.47037005]
        
        # # self.play(
        # #     angle_earth.animate.increment_value(2 * PI),
        # #     angle_moon.animate.increment_value(8 * 2 * PI),
        # #     run_time=20,
        # #     rate_func=linear
        # # )
        self.embed()
        

class CircleSurface(Surface):
    def uv_func(self, u, v):
        r = 4  # Raio máximo do círculo
        return np.array([
            r * np.cos(u) * v,
            r * np.sin(u) * v,
            0
        ])
    
class Teste2(Scene):
    def construct(self):
        self.camera.background_rgba = [0, 0, 0, 1]
        self.frame.reorient(43, 76, 1, IN, 10)
        sphere = Sphere(radius=2, depth_test=False)
        sphere.set_opacity(0.5)  # Visível mas transparente
        sphere.set_color(BLUE)
        surface = CircleSurface(
            u_range=[0, TAU],
            v_range=[0, 1],
            resolution=(100, 30),
            color=GREEN,
            depth_test=False
        )
        # surface.set_opacity(1)  # unnecesary, 1.0 is default
        self.add(sphere)
        self.play(ShowCreation(surface), run_time=3)
        self.embed()