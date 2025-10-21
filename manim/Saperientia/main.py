from manimlib import *
from manimlib.Saperientia.astro_structures import * 
from manimlib.Saperientia.stellarium import * 
from manimlib.Saperientia.Variables import * 
from manimlib.Saperientia.sap_rate_func import * 
from manimlib.Saperientia.camera import * 
from manimlib.Saperientia.video_mobject import *
import numpy as np
from scipy.spatial.transform import Rotation as R

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


#Exemplos  

class Demo(Scene):
    def construct(self):
        self.frame.reorient(0, 0, 0, ORIGIN+[0,0,0.001], 0.0005)
        # self.camera.frame.set_width(10000)
        stars,_=Stars()

        self.add(stars)
    
        sun = Sun()
        sun.add_updater(lambda m:self.add(sun))
        earth = Earth().move_to(RENDER_SUN_EARTH_DISTANCE*RIGHT)
        light = self.camera.light_source
        moon = Moon()  # Lua ainda menor

        self.add(sun, earth,moon)
        
    
        # Trackers para os ângulos da Terra e da Lua
        angle_earth = ValueTracker(0)
        angle_moon = ValueTracker(0)

        # Atualiza a posição da Terra orbitando o Sol
        earth_orbit = always_redraw(lambda: 
            earth.move_to(sun.get_center() + RENDER_SUN_EARTH_DISTANCE * np.array([
                np.cos(angle_earth.get_value()),
                np.sin(angle_earth.get_value()),
                0
            ]))
        )
        
        # Atualiza a posição da Lua orbitando a Terra
        moon_orbit = always_redraw(lambda:
            moon.move_to(earth.get_center() + RENDER_EARTH_MOON_DISTANCE * np.array([
                np.cos(angle_moon.get_value()),
                np.sin(angle_moon.get_value()),
                0
            ]))
        )
        # self.frame.reorient(330,70,0,[53.38888931,  0.        ,  0.        ],15)
        self.frame.reorient(336,70,0,earth.get_center(),1)
        esfera_celeste = EsferaCeleste().move_to(earth.get_center())
        equador = Equador(90).move_to(earth.get_center())
        eixopolar = EixoPolar(90).move_to(earth.get_center())
        esfera_celeste.add_updater(lambda m:self.add(esfera_celeste))
        self.play(ShowCreation(esfera_celeste),run_time=2)
        self.play(ShowCreation(equador),ShowCreation(eixopolar),run_time=4)
        self.play(Rotate(earth, 2*TAU,axis=Z_AXIS),run_time=4)
        p1 = PontoAstro(-30,200,center=earth.get_center(),cor=RED)
        p2 = PontoAstro(30,190,center=earth.get_center(),cor=GREEN)
        p3 = PontoAstro(10,240,center=earth.get_center(),cor=YELLOW)
        self.play(FadeIn(p1),FadeIn(p2),FadeIn(p3))
        arc1 = GrandeArco(p1,p2,cor=PINK, espessura=4,center=earth.get_center())
        arc2 = GrandeArco(p1,p3,cor=PINK, espessura=4,center=earth.get_center())
        arc3 = GrandeArco(p2,p3,cor=PINK, espessura=4,center=earth.get_center()) 
        
        self.play(ShowCreation(arc1),ShowCreation(arc2),ShowCreation(arc3))
        self.play(self.frame.animate.reorient(330,110,0,earth.get_center(),1),run_time=2)
        self.play(self.frame.animate.reorient(336,110,0,earth.get_center(),1),run_time=2)
        self.wait(2)
        self.play(self.frame.animate.reorient(330,70,0,earth.get_center(),1),run_time=2)
        ang1 = AnguloEsferico(p1,p3,p2,espessura=3,raio_circulo=0.7,center=earth.get_center())
        ang2 = AnguloEsferico(p3,p2,p1,espessura=3,raio_circulo=0.7,center=earth.get_center())
        ang3 = AnguloEsferico(p2,p1,p3,espessura=3,raio_circulo=0.7,center=earth.get_center())
        self.play(ShowCreation(ang1),ShowCreation(ang2),ShowCreation(ang3))
        self.play(
            angle_moon.animate.increment_value(4 * PI),
            self.frame.animate.reorient(60,80,0,earth.get_center(),3),
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

        self.embed()

    
        
