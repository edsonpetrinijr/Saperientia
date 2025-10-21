from manimlib import *
from manimlib.Saperientia.astro_structures import * 
from manimlib.Saperientia.stellarium import * 
from manimlib.Saperientia.Variables import * 
from manimlib.Saperientia.sap_rate_func import * 
from manimlib.Saperientia.camera import * 
from manimlib.Saperientia.video_mobject import *
import numpy as np
from scipy.spatial.transform import Rotation as R

#BODY 

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
        terra = Earth(radius=3)
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
        phi,theta,convert_camera_angles(80,180,0,90+23)
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
        
class Camera(Scene):
    def construct(self):
        earth = Earth(radius=RENDER_EARTH_RADIUS,clouds=False)
        self.add(earth)
        # earth.add_updater(lambda m,dt:m.rotate(dt*0.1,axis=OUT,about_point=ORIGIN))

        # Exemplo: quero “encostar” na origem, mantendo fovy do jeito que está
        frame = self.camera.frame

        # Exemplo: colar na origem, invertendo a direção
        reorientar_camera(
            frame,
            point=ORIGIN,
            theta_deg=0, phi_deg=90, gamma_deg=0,
        )
        self.wait(5)
        self.play(
            AnimarCamera(
                frame,
                point=earth.get_corner(OUT)*1.1,
                phi_deg=150,
                run_time=2
            )
        )

class NaSuper(Scene):
    def construct(self):
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
        
        earth = Earth(radius=RENDER_EARTH_RADIUS,clouds=False,resolution=(102,50)).rotate(180*DEGREES,axis=OUT,about_point=ORIGIN)
        self.add(earth)
       
        # Exemplo: quero “encostar” na origem, mantendo fovy do jeito que está
        frame = self.camera.frame

        # Exemplo: colar na origem, invertendo a direção

  
        p = P(0,0,raio=RENDER_EARTH_RADIUS)
        angle = ValueTracker(0)
        frame.add_updater(lambda m:reorientar_camera(m,np.array([np.sin(angle.get_value()),np.cos(angle.get_value()),0])*RENDER_EARTH_RADIUS*1.0000001,theta_deg=-angle.get_value()/PI*180,phi_deg=-20,fovy=70))
        self.play(angle.animate.increment_value(-2*PI),Rotate(earth,angle = 2*PI,axis=OUT,rate_func=linear,run_time=30),run_time=30,rate_func=linear)
        
# Helpers esféricos simples
def _latlon_to_xyz(lat_deg, lon_deg, R):
    lat = lat_deg * DEGREES
    lon = lon_deg * DEGREES
    x = R * np.cos(lat) * np.cos(lon)
    y = R * np.cos(lat) * np.sin(lon)
    z = R * np.sin(lat)
    return np.array([x, y, z], dtype=float)

class NaSuper2(Scene):
    def construct(self):
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

        earth = Earth(
            radius=RENDER_EARTH_RADIUS,
            clouds=False,
            resolution=(102, 50)
        ).rotate(180*DEGREES, axis=OUT, about_point=ORIGIN)
        self.add(earth)

        frame = self.camera.frame

        # Latitude alvo em graus. Exemplo São Paulo aproximado
        lat_deg = ValueTracker(-23.55)  
        R = RENDER_EARTH_RADIUS
        alt_factor = 1.0001  # levemente acima da superfície

        # Este ValueTracker governa a longitude. Use radianos como antes
        angle = ValueTracker(0)

        # Updater para qualquer latitude
        def camera_updater(mob):
            # longitude em graus a partir do seu angle em radianos
            lon_deg = angle.get_value() / PI * 180.0 * (-1.0)  # sinal negativo para casar com sua rotação
            # posição na superfície
            point = _latlon_to_xyz(lat_deg.get_value(), lon_deg, R * alt_factor)
            # ângulos pré flip olhando para o sul
            th_pre, ph_pre, ga_pre = lon_deg-90 , +lat_deg.get_value()-20 , 0
            # cola a câmera no ponto e ajusta orientação
            reorientar_camera(
                mob,
                point=point,
                theta_deg=th_pre,
                phi_deg=ph_pre,
                gamma_deg=ga_pre,
                epsilon=1e-6,
                fovy=60
            )

        frame.add_updater(camera_updater)

        # A Terra girando no mesmo sentido e velocidade angular coerente
        self.play(lat_deg.animate.increment_value(110),run_time=5)
        self.play(
            angle.animate.increment_value(-0.5*PI),
              # passa pelas longitudes conforme antes
            Rotate(earth, angle=0.5*PI, axis=OUT, rate_func=linear, run_time=10),
            run_time=10, rate_func=linear
        )
        frame.clear_updaters()
        self.play(frame.animate.reorient(0,70,0,earth.get_corner(OUT)*1.3,RENDER_EARTH_RADIUS*0.02),            
                  Rotate(earth, angle=0.5*PI, axis=OUT, rate_func=linear, run_time=2),run_time=2,rate_func=linear)
    
        self.play(frame.animate.reorient(0,70,0,earth.get_corner(OUT)*1.3,RENDER_EARTH_RADIUS*2),            
                  Rotate(earth, angle=0.5*PI, axis=OUT, rate_func=linear, run_time=4),run_time=4,rate_func=linear)
        lon_deg = angle.get_value() / PI * 180.0 * (-1.0)  # sinal negativo para casar com sua rotação
        # posição na superfície
        point = _latlon_to_xyz(lat_deg.get_value(), lon_deg, R * alt_factor)
        # ângulos pré flip olhando para o sul
        th_pre, ph_pre, ga_pre = lon_deg-90 , +lat_deg.get_value()-20 , 0
        # cola a câmera no ponto e ajusta orientação
    
        self.play(AnimarCamera(frame,point=point,theta_deg=th_pre,phi_deg=90,gamma_deg=ga_pre),run_time=4)
        self.play(AnimarCamera(frame,point=point,theta_deg=th_pre,phi_deg=ph_pre,gamma_deg=ga_pre),run_time=1)
        frame.add_updater(camera_updater)
        self.play(lat_deg.animate.increment_value(-70),run_time=5)
        self.play(
            angle.animate.increment_value(-0.5*PI),
              # passa pelas longitudes conforme antes
            Rotate(earth, angle=0.5*PI, axis=OUT, rate_func=linear, run_time=10),
            run_time=10, rate_func=linear
        )
        frame.clear_updaters()
        self.play(
            Rotate(earth, angle=2*PI, axis=OUT, rate_func=linear, run_time=10),
            AnimarCamera(frame,point=frame.get_center()*2,phi_deg=-120),
            run_time=10, rate_func=linear
        )
  

        self.wait(0.5)
        
class Observer(Scene):
    def construct(self):
        superficie=SuperficieObservador()
        self.add(superficie)
        equador = Equador(0,cor=BLUE_C,espessura=3*LINE_SIZE,)
        esfera = EsferaCeleste()
        self.add(esfera,equador)
        esfera.add_updater(lambda m: self.add(esfera))
        self.frame.reorient(90,86,0,ORIGIN+[0,0,0.1*RENDER_CELESTIAL_SPHERE_RADIUS],RENDER_CELESTIAL_SPHERE_RADIUS*0.1)
        direcao = PontosCardeais()
        self.wait(1)
        earth = Earth(radius=RENDER_CELESTIAL_SPHERE_RADIUS*10,clouds=False).rotate(-90*DEG,axis=RIGHT).align_to(superficie,OUT)
        equador2 = Equador(0,cor=BLUE_C,espessura=3*LINE_SIZE,raio=RENDER_CELESTIAL_SPHERE_RADIUS*100)
        equador3 = Equador(0,cor=BLUE_C,espessura=3*LINE_SIZE,raio=RENDER_CELESTIAL_SPHERE_RADIUS*10.1).move_to(earth.get_center())
        self.play(
            self.frame.animate.reorient(90,60,0,ORIGIN+[0,0,0.1*RENDER_CELESTIAL_SPHERE_RADIUS],RENDER_CELESTIAL_SPHERE_RADIUS*15),
            FadeIn(earth),
            FadeIn(equador2),FadeIn(equador3),run_time=3)
        group = Group()
        self.wait(4)
        eixo = EixoPolar(0,espessura=0.06*ELEMENTS_SCALE,comprimento=4*RENDER_CELESTIAL_SPHERE_RADIUS).set_z_index(20)
        eixo2 = EixoPolar(0,espessura=0.1*ELEMENTS_SCALE,comprimento=100*RENDER_CELESTIAL_SPHERE_RADIUS).set_z_index(2).move_to(earth.get_center())
        self.play(ShowCreation(eixo))
        self.play(ShowCreation(eixo2))


        group.add(superficie,esfera)
        eixo.add_updater(lambda m: self.add(eixo.move_to(group.get_center())))
        equador.add_updater(lambda m: m.move_to(eixo.get_center()))
        esfera.clear_updaters()
        esfera.add_updater(lambda m: self.add(esfera))
        esfera.deactivate_depth_test()
        self.play(
                  Rotate(group,15*DEG,RIGHT,about_point=earth.get_center()),run_time=10
                  )
        
        self.play(self.frame.animate.reorient(90,80,0,earth.get_center()+[0,0,4*RENDER_CELESTIAL_SPHERE_RADIUS],RENDER_CELESTIAL_SPHERE_RADIUS*20),run_time=2)
        self.wait(2)
        self.play(
            Rotate(group,-15*DEG,RIGHT,about_point=earth.get_center()),run_time=2
                )         
        self.play(
            Rotate(group,90*DEG,RIGHT,about_point=earth.get_center()),run_time=6
                )         
        self.play(self.frame.animate.reorient(30,80,0,earth.get_center()+[0,0,4*RENDER_CELESTIAL_SPHERE_RADIUS],RENDER_CELESTIAL_SPHERE_RADIUS*20),run_time=2)
#NESSE OBSERVER 2 DA PRA VER QUE A TEXTURA DA TERRA APARECE DO OUTRO LADO, TALVEZ DE PRA USAR ISSO DE TEXTURA PRA COLOCAR ESFERA CELESTE E TBM FAZER A ESFERA CELESTE 
        
class DevelopSky(Scene):
    def construct(self):
        stars,_ = Stars()
        self.add(stars)
        self.frame.reorient(0,70,0,ORIGIN,2)
        self.play(stars.animate.become(Stars(sphere_mode=True)[0]),run_time=7,rate_func=fast_to_slow)
        self.play(stars.animate.become(Stars(ayre=True)[0]),run_time=4)
        self.frame.add_ambient_rotation(0.5)
        self.wait(5)
        
class HRDiagram(Scene):
    def construct(self):
        self.frame.reorient(170,130,0,ORIGIN,3)
        axes = Axes().shift((LEFT+DOWN)*2)
        stars0,_ = Stars(hr_mode_radius=True,size_factor=4)
        stars1,_ = Stars(hr_mode=True,size_factor=4)
        stars2,_ = Stars(sphere_mode=True)
        stars3,_ = Stars()
        self.add(stars3)
        self.wait()
        self.play(self.frame.animate.reorient(0,40,0,ORIGIN,3),run_time=3)
        self.play(stars3.animate.become(stars2),run_time=7,rate_func=fast_to_slow)
        self.wait(3)
        self.play(stars3.animate.become(stars1),ShowCreation(axes),self.frame.animate.reorient(0,7,0,ORIGIN,6),run_time=3)
        self.wait(3)
        surface = ParametricSurface(
            lambda u, v: np.array([u,v, np.sqrt((10**((v+1.5)*2))*3.8e26/(4*PI*5.67e-8*(4600*((1/(0.92*((u/2)+0.656)+1.7))+(1/(0.92*((u/2)+0.656)+0.62))))**4))/695_500_000/15]),
            u_range=[-2, 6],
            v_range=[-2, 6],
            resolution=(50, 50),
        ).set_color(PURPLE).set_opacity(0.7).set_z_index(-3)
        self.play(stars3.animate.become(stars0),self.frame.animate.reorient(-30,90,0,ORIGIN+[0,0,2],8),run_time=3)
        self.wait(3)
        self.play(ShowCreation(surface))
        self.frame.add_ambient_rotation(0.5)
        self.wait(13)
        self.play(self.frame.animate.reorient(phi_degrees=150,gamma_degrees=0,center=ORIGIN+[0,0,2],height=6))
        self.wait(3)
        
class DevelopSky2(Scene):
    def construct(self):
        constellation = "UMA"
        self.look_at_constellation(constellation)
        self.frame.set_field_of_view(60*DEG)
        stars = Stars(time_years=-8000)
        stars2 = Stars(time_years=40000)
        stars3 = Stars()
        borders = Constellations(sphere_radius=10000,constellations=[constellation])
        borders_all = Constellations(sphere_radius=10000)
        asterisms = Asterisms(time_years=-8000,constellations=[constellation],sphere_radius=10000)
        asterisms_all = Asterisms(sphere_radius=10000)
        asterisms2 = Asterisms(time_years=40000,constellations=[constellation],sphere_radius=10000)
        asterisms3 = Asterisms(constellations=[constellation],sphere_radius=10000)
        self.add(stars)
        self.play(ShowCreation(borders),ShowCreation(asterisms),run_time=4)
        
        tracker = ValueTracker(-8000)
        text = Text(f"Ano {tracker.get_value():.0f}", font_size=48)
        text.to_corner(UL)
        text.fix_in_frame() 
        text.add_updater(
            lambda t: t.become(
                Text(f"Ano {tracker.get_value():.0f}", font_size=48)
                .to_corner(UL)
                .fix_in_frame()
            )
        )
        self.add(text)
        self.play(Write(text),run_time=2)
        self.play(stars.animate.become(stars2),asterisms.animate.become(asterisms2),tracker.animate.set_value(40000),run_time=10,rate_func=linear)
        self.wait(3)
        self.play(stars.animate.become(stars3),asterisms.animate.become(asterisms3),tracker.animate.set_value(2025),run_time=5,rate_func=linear)
        self.frame.add_ambient_rotation(0.5)
        self.play(self.frame.animate.reorient(0,90,0),FadeIn(borders_all),FadeIn(asterisms_all),FadeOut(text),FadeOut(asterisms),FadeOut(borders),run_time=3)
        self.wait(3)
        self.frame.clear_updaters()
        self.play(FadeOut(borders_all),FadeOut(asterisms_all),run_time=2)
        self.play(self.frame.animate.reorient(170,130,0,ORIGIN,3))
        axes = Axes().shift((LEFT+DOWN)*2)
        stars0 = Stars(hr_mode_radius=True,size_factor=4)
        stars1 = Stars(hr_mode=True,size_factor=4)
        stars2 = Stars(sphere_mode=True)
        self.wait()
        self.play(self.frame.animate.reorient(0,40,0,ORIGIN,3),run_time=3)
        self.play(stars.animate.become(stars2),run_time=7,rate_func=fast_to_slow)
        self.wait(3)
        self.play(stars.animate.become(stars1),ShowCreation(axes),self.frame.animate.reorient(0,7,0,ORIGIN,6),run_time=3)
        self.wait(3)
        surface = ParametricSurface(
            lambda u, v: np.array([u,v, np.sqrt((10**((v+1.5)*2))*3.8e26/(4*PI*5.67e-8*(4600*((1/(0.92*((u/2)+0.656)+1.7))+(1/(0.92*((u/2)+0.656)+0.62))))**4))/695_500_000/15]),
            u_range=[-2, 6],
            v_range=[-2, 6],
            resolution=(50, 50),
        ).set_color(PURPLE).set_opacity(0.7).set_z_index(-3)
        self.play(stars.animate.become(stars0),self.frame.animate.reorient(-30,90,0,ORIGIN+[0,0,2],8),run_time=3)
        self.wait(3)
        self.play(ShowCreation(surface))
        self.frame.add_ambient_rotation(0.5)
        self.wait(13)
        self.play(self.frame.animate.reorient(phi_degrees=150,gamma_degrees=0,center=ORIGIN+[0,0,2],height=6))
        self.wait(3)
        self.play(FadeOut(axes,rate_func=linear),FadeOut(surface,rate_func=linear),stars.animate.become(stars3),run_time=5,rate_func=slow_to_fast)
        self.wait(3)
        
        
class TesteEstrelas(Scene):
    def construct(self):
        self.frame.reorient(0,130,0,ORIGIN,1)
        estrelas = Stars()
        asterismos = Asterisms(sphere_radius=100000)
        constela = Constellations(sphere_radius=100000)
        
        self.add(estrelas,constela,asterismos)
