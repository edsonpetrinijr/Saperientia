from manimlib import *
from manimlib.Saperientia.astro_structures import * 
from manimlib.Saperientia.stellarium import * 
from manimlib.Saperientia.Variables import * 
from manimlib.Saperientia.sap_rate_func import * 
from manimlib.Saperientia.camera import * 
from manimlib.Saperientia.video_mobject import *
import numpy as np
from scipy.spatial.transform import Rotation as R



class DiaSideral(Scene):
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
        self.frame.reorient(90,80,0,ORIGIN,CELESTIAL_SPHERE_RADIUS*1)
        terra = Earth(clouds=False)
        self.add(terra)
        terra.add_updater(lambda m,dt: m.rotate(0.6*dt,axis=OUT))
        self.wait(15)
        self.play(self.frame.animate.reorient(90,0,0,ORIGIN,CELESTIAL_SPHERE_RADIUS*1),run_time=2)
        linha = Line(start=[0,0,0],end=terra.get_right()*4)
        self.play(ShowCreation(linha))
        linha.add_updater(lambda m,dt:m.rotate(0.6*dt,axis=OUT,about_point=ORIGIN))
        self.wait(12)
        self.embed()
        
class Pendulo(Scene):
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
        self.frame.reorient(90,80,0,ORIGIN,CELESTIAL_SPHERE_RADIUS*1)
        terra = Earth(clouds=False)
        self.add(terra)
        terra.add_updater(lambda m,dt: m.rotate(0.6*dt,axis=OUT))
        L = 0.01
        g = 0.1
        theta0 = 1
        gamma = 0
        omega = np.sqrt(g / L)
        bob_radius = 0.001

        # ----- Função de interpolação (oscilador amortecido) -----
        def theta_of_time(t):
            return theta0 * np.cos(omega * t) * np.exp(-gamma * t)

        # ----- Elementos da cena -----
        # Agora consideramos Z como eixo vertical (para cima)
        pivot = Sphere(radius=0.0005, color=PURPLE).shift(OUT * 0)  # pivô em (0,0,0)
        bob = Sphere(radius=bob_radius)

        # O fio inicial vai do pivô até L unidades abaixo (eixo Z negativo)
        rod = Line(ORIGIN, -L * OUT, stroke_width=3)

        # suporte no topo

        pendulum_group = Group(pivot, rod, bob)

        # posição inicial (movimento no plano XZ)
        initial_theta = theta_of_time(0)
        bob_pos = np.array([
            L * np.sin(initial_theta),  # x
            0,                          # y (plano da tela)
            -L * np.cos(initial_theta)  # z para baixo
        ])
        rod.put_start_and_end_on(pivot.get_center(), bob_pos)
        bob.move_to(bob_pos)

        self.add(pendulum_group)
        des=RENDER_EARTH_RADIUS*OUT*1.2
        pendulum_group.shift(des)
        # ----- Tracker de tempo -----
        t_tracker = ValueTracker(0.0)

        def advance_time(mob, dt):
            mob.set_value(mob.get_value() + dt)
        t_tracker.add_updater(advance_time)

        # ----- Updater de movimento -----
        def update_pendulum(mob):
            t = t_tracker.get_value()
            th = theta_of_time(t)
            pos = np.array([
                L * np.sin(th),
                0,
                -L * np.cos(th)
            ]+des)
            rod.put_start_and_end_on(pivot.get_center(), pos)
            bob.move_to(pos)

        rod.add_updater(lambda m: update_pendulum(m))
        bob.add_updater(lambda m: update_pendulum(m))

        self.add(t_tracker)

        # ----- Câmera (olhando de cima, com Z para cima) -----
        self.wait(5)
        self.play(self.frame.animate.reorient(180,60,0,terra.get_corner(OUT),CELESTIAL_SPHERE_RADIUS*0.3))
        self.wait(10)
        self.play(self.frame.animate.reorient(180,0,0,terra.get_corner(OUT),CELESTIAL_SPHERE_RADIUS*0.3))
        self.wait(10)

        # Limpeza
        t_tracker.remove_updater(advance_time)
        rod.clear_updaters()
        bob.clear_updaters()
        self.wait(2)
        self.embed()
        
class Eixo(Scene):
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
        self.frame.reorient(90,90,0,ORIGIN,CELESTIAL_SPHERE_RADIUS*1)
        terra = Earth(clouds=False)
        self.add(terra)
        terra.add_updater(lambda m,dt: m.rotate(0.6*dt,axis=OUT))
        self.wait(15)
        eixo = EixoPolar(90,espessura=0.02*ELEMENTS_SCALE)
        p1 = P(90,0,cor=PURPLE,raio=RENDER_EARTH_RADIUS).scale(0.4)
        p2 = P(-90,0,cor=PURPLE,raio=RENDER_EARTH_RADIUS).scale(0.4)
        self.play(ShowCreation(p1),ShowCreation(p2))
        self.play(self.frame.animate.reorient(90,10,0,ORIGIN,CELESTIAL_SPHERE_RADIUS*0.5),run_time=2)
        self.wait(5)
        self.play(self.frame.animate.reorient(100,170,0,ORIGIN,CELESTIAL_SPHERE_RADIUS*0.5),run_time=2)
        self.wait(5)
        self.play(self.frame.animate.reorient(90,90,0,ORIGIN,CELESTIAL_SPHERE_RADIUS*2),run_time=2)
        self.play(ShowCreation(eixo))
        self.wait(5)
        label1 = Text("Norte",font_size=3).scale(1).rotate(90*DEGREES,axis=RIGHT).rotate(90*DEGREES,axis=OUT).next_to(p1.get_center(),OUT,buff=0.01)
        label2 = Text("Sul",font_size=3).scale(1).rotate(90*DEGREES,axis=RIGHT).rotate(90*DEGREES,axis=OUT).next_to(p2.get_center(),IN,buff=0.01)
        self.play(Write(label1))
        self.play(Write(label2))
        self.wait(10)
        terra.clear_updaters()
        self.play(terra.animate.become(Earth(clouds=False).rotate(165*DEGREES,axis=OUT)))
        self.play(self.frame.animate.reorient(120,90,0,ORIGIN,CELESTIAL_SPHERE_RADIUS*1),run_time=2)
        self.play(FadeOut(label1),FadeOut(label2))
        pessoa = PontoAstro(0,90,cor=BLUE,raio=RENDER_EARTH_RADIUS).scale(0.4)
        self.play(FadeIn(pessoa))
        self.wait(2)
        arco = GrandeArco(pessoa,p1,cor=GREEN,espessura=LINE_SIZE*2)
        self.play(ShowCreation(arco))
        seta= Seta(pessoa,0,cor=GREEN,espessura=LINE_SIZE*2,comprimento=0.01)
        self.wait(5)
        self.play(ShowCreation(seta))
        self.wait(5)
        self.play(self.frame.animate.reorient(120,20,0,ORIGIN,CELESTIAL_SPHERE_RADIUS*1),run_time=2)
        self.wait(5)
        self.play(self.frame.animate.reorient(120,90,0,ORIGIN,CELESTIAL_SPHERE_RADIUS*1),run_time=2)
        self.wait(2)
        grade = Grade(raio=RENDER_EARTH_RADIUS,cor_ar=GREEN).set_z_index(10)
        self.play(ShowCreation(grade[1]))
        self.play(terra.animate.set_opacity(0.6))
        
class SolarSystem(Scene):
    def construct(self):
        stars_data = extract_star_data()
        stars = Group()
        for x, y, z, size, color in stars_data:
            star = DotCloud(
                points=[2*np.array([x, y, z])],
                color=color,
                radius=size,
                opacity=0.75
            )
            stars.add(star)
        self.add(stars)
        
        sun=Sun(radius=RENDER_SUN_RADIUS*5)
        mercury = Mercury().shift(RIGHT * RENDER_SUN_MERCURY_DISTANCE)
        venus = Venus().shift(RIGHT * RENDER_SUN_VENUS_DISTANCE)
        earth = Earth().shift(RIGHT * RENDER_SUN_EARTH_DISTANCE)
        mars = Mars().shift(RIGHT * RENDER_SUN_MARS_DISTANCE)
        jupiter = Jupiter().shift(RIGHT * RENDER_SUN_JUPITER_DISTANCE)
        saturn = Saturn().shift(RIGHT * RENDER_SUN_SATURN_DISTANCE)
        uranus = Uranus().shift(RIGHT * RENDER_SUN_URANUS_DISTANCE)
        neptune = Neptune().shift(RIGHT * RENDER_SUN_NEPTUNE_DISTANCE)

        # ========== OPCIONAL: ÓRBITAS ==========
        orbits = [
            Circle(radius=RENDER_SUN_MERCURY_DISTANCE, color=GREY_A, stroke_opacity=0.2),
            Circle(radius=RENDER_SUN_VENUS_DISTANCE, color=GREY_A, stroke_opacity=0.2),
            Circle(radius=RENDER_SUN_EARTH_DISTANCE, color=GREY_A, stroke_opacity=0.2),
            Circle(radius=RENDER_SUN_MARS_DISTANCE, color=GREY_A, stroke_opacity=0.2),
            Circle(radius=RENDER_SUN_JUPITER_DISTANCE, color=GREY_A, stroke_opacity=0.2),
            Circle(radius=RENDER_SUN_SATURN_DISTANCE, color=GREY_A, stroke_opacity=0.2),
            Circle(radius=RENDER_SUN_URANUS_DISTANCE, color=GREY_A, stroke_opacity=0.2),
            Circle(radius=RENDER_SUN_NEPTUNE_DISTANCE, color=GREY_A, stroke_opacity=0.2),
        ]
        for orbit in orbits:
            self.add(orbit)

        # ========== ADICIONA OS CORPOS NA CENA ==========
        self.add(
            sun,
            mercury, venus, earth, mars,
            jupiter, saturn, uranus, neptune
        )
        eixo=EixoPolar(90).move_to(earth.get_center())
        self.add(eixo)
        earth.add_updater(lambda m,dt:m.rotate(dt*0.1,axis=eixo.get_start()-eixo.get_end(),about_point=m.get_center()))
        terra = Group(eixo,earth)
        terra.rotate(OBLIQUIDADE*DEG,axis=UP)
        self.frame.reorient(0,80,0,earth.get_center(),CELESTIAL_SPHERE_RADIUS*3),
        self.wait(5)
        self.play(self.frame.animate.reorient(0,80,0,earth.get_center(),CELESTIAL_SPHERE_RADIUS*30),run_time=2)
        self.wait()
        self.play(self.frame.animate.reorient(0,83,0,(saturn.get_center()+uranus.get_center())/2,CELESTIAL_SPHERE_RADIUS*1300),run_time=2)
        self.play(eixo.animate.scale(100))
        
        obliquities = {
            mercury: MERCURY_OBLIQUITY,
            venus: VENUS_OBLIQUITY,
            earth: EARTH_OBLIQUITY,
            mars: MARS_OBLIQUITY,
            jupiter: JUPITER_OBLIQUITY,
            saturn: SATURN_OBLIQUITY,
            uranus: URANUS_OBLIQUITY,
            neptune: NEPTUNE_OBLIQUITY,
        }
        # ===== CRIA EIXOS E ROTACIONA =====
        planet_groups = []
        eixos = Group()
        for planet, obliq in obliquities.items():
            eixo = EixoPolar(90,cor=BLUE).move_to(planet.get_center())
            self.add(eixo)
            planet.rotate(obliq * DEG, axis=UP)
            eixo.rotate(obliq * DEG, axis=UP)
            eixos.add(eixo.copy())
            
        eixo1=eixos[0]
        eixo2=eixos[1]
        eixo3=eixos[2]
        eixo4=eixos[3]
        eixo5=eixos[4]
        eixo6=eixos[5]
        eixo7=eixos[6]
        eixo8=eixos[7]
        
        self.play(eixo1.animate.scale(100),eixo2.animate.scale(100),eixo3.animate.scale(100),eixo4.animate.scale(100),eixo5.animate.scale(100),eixo6.animate.scale(100),
                  eixo7.animate.scale(100),eixo8.animate.scale(100),run_time=2
                  )
        self.wait(3)
        plano = PlanoRetangular(3*RENDER_SUN_NEPTUNE_DISTANCE,3*RENDER_SUN_NEPTUNE_DISTANCE,cor=PURPLE).set_opacity(0.3)
        self.play(ShowCreation(plano))
        self.play(self.frame.animate.reorient(0,90,0,(saturn.get_center()+uranus.get_center())/2,CELESTIAL_SPHERE_RADIUS*1300))
        self.wait(3)
        self.camera.light_source.move_to(mercury.get_center())
        self.play(self.frame.animate.reorient(0,85,0,venus.get_center(),CELESTIAL_SPHERE_RADIUS*1.7),eixo2.animate.scale(0.01),run_time=2)
        venus.add_updater(lambda m,dt:m.rotate(dt*0.2,axis=eixo2.get_start()-eixo2.get_end(),about_point=m.get_center()))
        self.wait(5)
        self.frame.add_ambient_rotation(1*DEG)
        self.play(self.frame.animate.reorient(180,83,30,(saturn.get_center()+uranus.get_center())/2,CELESTIAL_SPHERE_RADIUS*1300),run_time=6)
        
class LesteOeste(Scene):
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
        self.frame.reorient(90,85,0,ORIGIN,CELESTIAL_SPHERE_RADIUS*1)
        terra = Earth(clouds=False)
        
        self.add(terra)
        pessoa = PontoAstro(0,90,cor=BLUE,raio=RENDER_EARTH_RADIUS).scale(0.4)
        self.play(FadeIn(pessoa))
        terra.add_updater(lambda m,dt: m.rotate(0.6*dt,axis=OUT))
        pessoa.add_updater(lambda m,dt: m.rotate(0.6*dt,axis=OUT,about_point=ORIGIN))
        self.wait(10)
        
        seta_leste = Arrow(ORIGIN,2*RIGHT).fix_in_frame()
        seta_oeste = Arrow(ORIGIN,-2*RIGHT).fix_in_frame()
        seta_norte = Arrow(ORIGIN,2*UP).fix_in_frame()
        seta_sul = Arrow(ORIGIN,-2*UP).fix_in_frame()
        oeste=Text('Oeste').shift(LEFT*5).fix_in_frame()
        leste=Text('Leste').shift(RIGHT*5).fix_in_frame()
        self.play(FadeIn(seta_leste))
        self.wait(8)
        terra.clear_updaters()
        pessoa.clear_updaters()
        self.play(
            FadeOut(seta_leste),
            FadeOut(pessoa),
            run_time=2
        )
        pessoa2 = PontoAstro(0,270,cor=BLUE,raio=RENDER_EARTH_RADIUS).scale(0.4)
        pessoa3 = PontoAstro(30,270,cor=BLUE,raio=RENDER_EARTH_RADIUS).scale(0.4)
        pessoa4 = PontoAstro(-30,270,cor=BLUE,raio=RENDER_EARTH_RADIUS).scale(0.4)
        pessoa5 = PontoAstro(60,270,cor=BLUE,raio=RENDER_EARTH_RADIUS).scale(0.4)
        pessoa6 = PontoAstro(-60,270,cor=BLUE,raio=RENDER_EARTH_RADIUS).scale(0.4)
        self.play(
            FadeIn(pessoa2),
            self.frame.animate.reorient(270,85,0,pessoa2.get_center(),CELESTIAL_SPHERE_RADIUS*0.7))
        semi=GrandeCirculo([0,1,0],raio=RENDER_EARTH_RADIUS*1.001,espessura=1.7*LINE_SIZE)
        self.play(ShowCreation(semi))
        self.play(
            FadeIn(pessoa3),
            FadeIn(pessoa4),
            FadeIn(pessoa5),
            FadeIn(pessoa6),
        )
        setal1= Seta(pessoa2,-90,cor=GREEN,espessura=LINE_SIZE*1.5,comprimento=0.02)
        setal2= Seta(pessoa3,-90,cor=GREEN,espessura=LINE_SIZE*1.5,comprimento=0.02)
        setal3= Seta(pessoa4,-90,cor=GREEN,espessura=LINE_SIZE*1.5,comprimento=0.02)
        setal4= Seta(pessoa5,-90,cor=GREEN,espessura=LINE_SIZE*1.5,comprimento=0.02)
        setal5= Seta(pessoa6,-90,cor=GREEN,espessura=LINE_SIZE*1.5,comprimento=0.02)
        self.play(ShowCreation(setal1),ShowCreation(setal2),ShowCreation(setal3),ShowCreation(setal4),ShowCreation(setal5))
        self.wait()
        setao1= Seta(pessoa2,90,cor=GREEN,espessura=LINE_SIZE*1.5,comprimento=0.02)
        setao2= Seta(pessoa3,90,cor=GREEN,espessura=LINE_SIZE*1.5,comprimento=0.02)
        setao3= Seta(pessoa4,90,cor=GREEN,espessura=LINE_SIZE*1.5,comprimento=0.02)
        setao4= Seta(pessoa5,90,cor=GREEN,espessura=LINE_SIZE*1.5,comprimento=0.02)
        setao5= Seta(pessoa6,90,cor=GREEN,espessura=LINE_SIZE*1.5,comprimento=0.02)
        self.play(ShowCreation(setao1),ShowCreation(setao2),ShowCreation(setao3),ShowCreation(setao4),ShowCreation(setao5))
        self.wait(1)
        self.play(ShowCreation(seta_leste))
        self.play(ShowCreation(seta_oeste))
        self.play(ShowCreation(seta_norte))
        self.play(ShowCreation(seta_sul))
        self.wait(1)
        self.play(self.frame.animate.reorient(270,60,0,pessoa5.get_center(),CELESTIAL_SPHERE_RADIUS*0.7),run_time=3)
        self.wait(2)
        self.play(self.frame.animate.reorient(270,150,0,pessoa6.get_center(),CELESTIAL_SPHERE_RADIUS*0.7),run_time=3)
        self.wait(2)
        self.play(self.frame.animate.reorient(270,85,0,pessoa2.get_center(),CELESTIAL_SPHERE_RADIUS*0.7))
        self.wait(5)
        self.play(
            FadeOut(seta_leste),
            FadeOut(seta_oeste),
            FadeOut(seta_norte),
            FadeOut(seta_sul)
        )
        self.play(terra.animate.set_opacity(0.6))
        grade = Grade(raio=RENDER_EARTH_RADIUS,cor_dec=BLUE).set_z_index(10)
        self.play(ShowCreation(grade[0]))
        self.wait()
        self.play(self.frame.animate.reorient(270,85,0,ORIGIN,CELESTIAL_SPHERE_RADIUS*1.3))
        self.frame.add_ambient_rotation(2)
        self.wait(10)
        self.embed()

class LesteOesteMovimento(Scene):
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
        self.frame.reorient(90,50,0,ORIGIN,CELESTIAL_SPHERE_RADIUS*1)
        terra = Earth(clouds=False)
        
        self.add(terra)
        pessoa = PontoAstro(0,-120,cor=BLUE,raio=RENDER_EARTH_RADIUS).scale(0.4)
        self.play(FadeIn(pessoa))
        terra.add_updater(lambda m,dt: m.rotate(0.2*dt,axis=OUT))
        pessoa.add_updater(lambda m,dt: m.rotate(0.2*dt,axis=OUT,about_point=ORIGIN))
        self.wait(10)
        self.play(self.frame.animate.add_ambient_rotation(0.2),run_time=4)
        self.frame.add_ambient_rotation(0.2)
        self.wait(10)
        
class LesteOesteMovimento2(Scene):
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
        self.frame.reorient(270,90,0,ORIGIN,CELESTIAL_SPHERE_RADIUS*1)
        terra = Earth(clouds=False)
        sun = Sun().shift(LEFT*100)
        self.camera.light_source.move_to(LEFT*100)
        self.add(sun)
        
        self.add(terra)
        pessoa = PontoAstro(0,-120,cor=BLUE,raio=RENDER_EARTH_RADIUS).scale(0.4)
        self.play(FadeIn(pessoa))
        terra.add_updater(lambda m,dt: m.rotate(0.2*dt,axis=OUT))
        pessoa.add_updater(lambda m,dt: m.rotate(0.2*dt,axis=OUT,about_point=ORIGIN))
        self.wait(2)
        self.play(self.frame.animate.add_ambient_rotation(0.2),run_time=4)
        self.frame.add_ambient_rotation(0.2)
        self.wait(5)
        self.play(self.frame.animate.reorient(phi_degrees=90,gamma_degrees=-90,center=[5.9182663 , 2.37760878 ,0],height=CELESTIAL_SPHERE_RADIUS*0.7),run_time=4)
        self.frame.add_updater(lambda m:m.reorient(phi_degrees=90,gamma_degrees=-90,center=pessoa.get_center(),height=CELESTIAL_SPHERE_RADIUS*0.7))

        self.wait(15)

class Achatamento(Scene):
    def construct(self):
        self.frame.reorient(270,90,0,ORIGIN,CELESTIAL_SPHERE_RADIUS*1)
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
        earth = Earth()
        earth.add_updater(lambda m,dt: m.rotate(0.4*dt,axis=OUT))
        self.wait(4)
        self.play(earth.animate.apply_matrix(np.diag([1.1,1.1,0.9])),run_time=4)
        self.wait(4)
        equador = Equador(90,raio=RENDER_EARTH_RADIUS*1.15,espessura=6)
        self.play(FadeIn(equador))
        self.wait(1)
        
class EquadorDemonstra(Scene):
    def construct(self):
        self.frame.reorient(270,75,0,ORIGIN,CELESTIAL_SPHERE_RADIUS*1)
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
        earth = Earth(clouds=False)
        self.add(earth)
        self.play(earth.animate.set_opacity(0.5))
        earth.add_updater(lambda m:self.add(earth))
        pessoa = P(0,0).move_to(ORIGIN)
        self.play(FadeIn(pessoa))
        
        lugar_geometrico_1 = GrandeCirculo([0,0,1],raio=RENDER_EARTH_RADIUS,espessura=5)
        angle = ValueTracker(0)
    
        #Qualquer direação 
        marcador =  MarcadorAngulo(P(90,0,raio=RENDER_EARTH_RADIUS),lugar_geometrico_1.point_from_proportion(angle.get_value()/360),raio_arco=0.05)
        marcador.add_updater(lambda m: m.become(MarcadorAngulo(P(90,0,raio=RENDER_EARTH_RADIUS),lugar_geometrico_1.point_from_proportion(angle.get_value()/360),raio_arco=0.05)))
        self.play(angle.animate.increment_value(360),run_time=5,rate_func=linear)
        
        #Marcador
        self.add(marcador)
        self.frame.add_ambient_rotation(-0.7)
        #equador
        self.play(ShowCreation(lugar_geometrico_1),angle.animate.increment_value(360),run_time=8,rate_func=linear)
        self.wait(5)
        
class Temporeal(Scene):
    def construct(self):
        self.frame.reorient(270,90,0,ORIGIN,CELESTIAL_SPHERE_RADIUS*1)
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
        earth = Earth(clouds=False)
        self.play(earth.animate.rotate(0.25*DEG,axis=OUT),run_time=60)

class Circumpolar(Scene):
    def construct(self):
        self.frame.reorient(0,90,0,ORIGIN,CELESTIAL_SPHERE_RADIUS*1.5)
        stars_data = extract_star_data(send_name=True)
        stars = Group()
        star_lookup = {}
        for star_dic in stars_data:
            x, y, z = star_dic["x"], star_dic["y"], star_dic["z"]
            size = star_dic["size"]
            color = star_dic["color"]
            hip = star_dic["hip"]
            distance = np.linalg.norm([x, y, z])
            
            star_lookup[hip] = star_dic
            star = GlowDots(
                points=[np.array([x, y, z])*60/distance],
                color=color,
                radius=size*200/distance,
                opacity=0.75
            )
            
            stars.add(star)
        self.add(stars)
        ASTERISM_LINES = [
            (11767, 77055),  #Ursa Menor
            (79822, 77055),  #Ursa Menor
            (72607, 77055),  #Ursa Menor
            (75097, 72607),  #Ursa Menor
            (75097, 79822),  #Ursa Menor
            (67301, 65378),  # Alkaid – Mizar
            (65378, 62956),  # Mizar – Alioth
            (62956, 59774),  # Alioth – Megrez
            (59774, 54061),  # Megrez – Dubhe
            (59774, 58001),  # Megrez – Phecda
            (58001, 53910),  # Phecda – Merak
            (53910, 54061),  # Merak – Dubhe
            (60718, 61084), #Crux
            (62434, 59747), # Crux
            (107089,112405), #Oct
            (70638,107089), #Oct
            (70638,112405) #oct
            
  


        ]

        # === 3. Criar as linhas dos asterismos ===
        lines = Group()
        for hip1, hip2 in ASTERISM_LINES:
            if hip1 not in star_lookup or hip2 not in star_lookup:
                continue

            s1 = star_lookup[hip1]
            s2 = star_lookup[hip2]

            p1 = np.array([s1["x"], s1["y"], s1["z"]]) * 60 / np.linalg.norm([s1["x"], s1["y"], s1["z"]])
            p2 = np.array([s2["x"], s2["y"], s2["z"]]) * 60 / np.linalg.norm([s2["x"], s2["y"], s2["z"]])

            line = Line3D(
                p1, p2,
                color=WHITE,
                opacity=0.7
            )
            lines.add(line)


        earth = Earth(clouds=False)
        self.add(earth)
        eixo = EixoPolar(90,comprimento=140,espessura=0.04*ELEMENTS_SCALE)

        earth.add_updater(lambda m,dt: m.rotate(0.1*dt,axis=OUT))
        self.wait(2)
        self.play(self.frame.animate.reorient(0,156,0,ORIGIN,CELESTIAL_SPHERE_RADIUS*1.5),run_time=2)
        self.play(ShowCreation(eixo))
        self.wait(2)
        self.play(self.frame.animate.add_ambient_rotation(0.1),run_time=2)
        self.frame.add_ambient_rotation(0.1)
        self.wait(2)
        earth.add_updater(lambda m:self.add(earth))
        self.play(earth.animate.set_opacity(0.4))
        self.play(self.frame.animate.reorient(phi_degrees=175,gamma_degrees=0,center=ORIGIN+[0,0,16],height=CELESTIAL_SPHERE_RADIUS*1.2),run_time=4)
        self.wait(6)
        self.play(self.frame.animate.reorient(phi_degrees=90,gamma_degrees=0,center=ORIGIN,height=CELESTIAL_SPHERE_RADIUS*8),run_time=4)
        self.wait(28)
        self.play(earth.animate.set_opacity(1))
        
        self.play(self.frame.animate.reorient(phi_degrees=170,gamma_degrees=0,center=ORIGIN+[0,0,58],height=56),run_time=4)
        self.wait(5)
        self.play(ShowCreation(lines[:5]))
        self.wait(3)
        self.play(ShowCreation(lines[5:12]))
        self.wait(4)
        pessoa = P(0,270,raio=RENDER_EARTH_RADIUS)
        pessoa.add_updater(lambda m,dt:m.rotate(-0.1*dt,axis=OUT,about_point=ORIGIN))
        self.frame.add_updater(lambda m:m.reorient(phi_degrees=90,gamma_degrees=90,center=pessoa.get_center(),height=5))
        self.wait(5)
        self.frame.clear_updaters()
        self.frame.add_ambient_rotation(0.1)
        self.play(self.frame.animate.reorient(phi_degrees=90,gamma_degrees=90,center=pessoa.get_center(),height=65),run_time=4)
        self.wait(3)
        self.play(self.frame.animate.reorient(phi_degrees=10,gamma_degrees=180,center=ORIGIN+[0,0,-58],height=56),run_time=4)
        self.wait(12)
        self.play(ShowCreation(lines[12:14]))
        self.frame.clear_updaters()
        earth.clear_updaters()
        self.wait(5)
        dec = -87.31
        ra = 93.75
        theta = np.pi / 2 - np.radians(dec)
        phi = np.radians(ra) - np.pi / 2
        x=np.sin(theta)*np.sin(phi)*60
        y=np.sin(theta)*np.cos(phi)*60
        z=np.cos(theta)*60
        p1 = P(0,0).move_to(lines[12].get_start())
        p2 = P(0,0).move_to([x,y,z])
        arco = GrandeArco(p1,p2,espessura=4).set_opacity(0.7)
        self.play(ShowCreation(arco))
        self.wait(5)
        self.play(self.frame.animate.reorient(phi_degrees=10,gamma_degrees=180,center=ORIGIN+[1,1,-58],height=4),run_time=4)
        self.wait(3)
        self.play(self.frame.animate.reorient(phi_degrees=10,gamma_degrees=180,center=ORIGIN+[0,0,-58],height=56),run_time=4)
        self.play(ShowCreation(lines[14:]))
        
        
        self.wait(2)

class Precessao(Scene):
    def construct(self):
        self.frame.reorient(90,45,0,ORIGIN,CELESTIAL_SPHERE_RADIUS*3)
        light = self.camera.light_source
        terra = Earth(clouds=False)
        esfera_celeste = EsferaCeleste()
        self.add(terra,esfera_celeste)
        esfera_celeste.add_updater(lambda m: self.add(esfera_celeste))
        
        eixo = EixoPolar(90)
        equador = Equador(90,cor=YELLOW)
        self.add(eixo,equador)
        group = Group(terra,eixo,equador)
        self.wait(3)
        self.play(group.animate.rotate(23.5*DEGREES,axis=RIGHT,about_point=ORIGIN))
        ecliptica = Equador(90,cor=PINK)
        paralelo = Paralelo((90-23.5),90)
        self.play(ShowCreation(paralelo))
        terra.add_updater(lambda m,dt: m.rotate(-dt*0.3,axis=OUT,about_point=ORIGIN))
        eixo.add_updater(lambda m,dt: m.rotate(-dt*0.3,axis=OUT,about_point=ORIGIN))
        equador.add_updater(lambda m,dt: m.rotate(-dt*0.3,axis=OUT,about_point=ORIGIN))
        terra.add_updater(lambda m,dt: m.rotate(dt*20,axis=(eixo.get_end()-eixo.get_start())))
        self.wait(10)
        self.embed()   
class PrecessaoNutacao(Scene):
    def construct(self):
        self.frame.reorient(90,45,0,ORIGIN,CELESTIAL_SPHERE_RADIUS*3)
        light = self.camera.light_source
        terra = Earth(clouds=False).rotate(20*DEG,axis=LEFT)
        esfera_celeste = EsferaCeleste()
        self.add(terra,esfera_celeste)
        esfera_celeste.add_updater(lambda m: self.add(esfera_celeste))
        
        eixo = EixoPolar(90)
        equador = Equador(90,cor=YELLOW)
        self.add(eixo,equador)
        self.wait(3)
        ecliptica = Equador(90,cor=PINK)
        paralelo = Paralelo((90-23.5),90)
        self.play(ShowCreation(paralelo))
        def func_prec(t):
            p =P(90-23.5-np.sin(t/180*20*PI)*3,t)
            return p.get_center()
        t= ValueTracker(0)
        terra.add_updater(lambda m,dt: m.rotate(-dt*0.3,axis=OUT,about_point=ORIGIN))
        eixo.add_updater(lambda m: m.put_start_and_end_on(func_prec(t.get_value()),-func_prec(t.get_value())))
        equador.add_updater(lambda m,dt: m.become(GrandeCirculo(func_prec(t.get_value()))))
        terra.add_updater(lambda m,dt: m.rotate(dt*20,axis=(eixo.get_end()-eixo.get_start())))
        self.play(t.animate.increment_value(360),rate_func=linear,run_time=10)
        self.embed()   
        

class SolarSystemAxis(Scene):
    def construct(self):
        stars_data = extract_star_data()
        stars = Group()
        for x, y, z, size, color in stars_data:
            star = DotCloud(
                points=[2*np.array([x, y, z])],
                color=color,
                radius=size,
                opacity=0.75
            )
            stars.add(star)
        self.add(stars)
        
        sun=Sun(radius=RENDER_SUN_RADIUS*5)
        mercury = Mercury().shift(RIGHT * RENDER_SUN_MERCURY_DISTANCE)
        venus = Venus().shift(RIGHT * RENDER_SUN_VENUS_DISTANCE)
        earth = Earth(clouds=False).shift(RIGHT * RENDER_SUN_EARTH_DISTANCE)
        mars = Mars().shift(RIGHT * RENDER_SUN_MARS_DISTANCE)
        jupiter = Jupiter().shift(RIGHT * RENDER_SUN_JUPITER_DISTANCE)
        saturn = Saturn().shift(RIGHT * RENDER_SUN_SATURN_DISTANCE)
        uranus = Uranus().shift(RIGHT * RENDER_SUN_URANUS_DISTANCE)
        neptune = Neptune().shift(RIGHT * RENDER_SUN_NEPTUNE_DISTANCE)

        # ========== OPCIONAL: ÓRBITAS ==========
        orbits = [
            Circle(radius=RENDER_SUN_MERCURY_DISTANCE, color=GREY_A, stroke_opacity=0.2),
            Circle(radius=RENDER_SUN_VENUS_DISTANCE, color=GREY_A, stroke_opacity=0.2),
            Circle(radius=RENDER_SUN_EARTH_DISTANCE, color=GREY_A, stroke_opacity=0.2),
            Circle(radius=RENDER_SUN_MARS_DISTANCE, color=GREY_A, stroke_opacity=0.2),
            Circle(radius=RENDER_SUN_JUPITER_DISTANCE, color=GREY_A, stroke_opacity=0.2),
            Circle(radius=RENDER_SUN_SATURN_DISTANCE, color=GREY_A, stroke_opacity=0.2),
            Circle(radius=RENDER_SUN_URANUS_DISTANCE, color=GREY_A, stroke_opacity=0.2),
            Circle(radius=RENDER_SUN_NEPTUNE_DISTANCE, color=GREY_A, stroke_opacity=0.2),
        ]
        for orbit in orbits:
            self.add(orbit)

        # ========== ADICIONA OS CORPOS NA CENA ==========
        self.add(
            sun,
            mercury, venus, earth, mars,
            jupiter, saturn, uranus, neptune
        )
        eixo=EixoPolar(90).move_to(earth.get_center())
        self.add(eixo)
        earth.add_updater(lambda m,dt:m.rotate(dt*0.1,axis=eixo.get_start()-eixo.get_end(),about_point=m.get_center()))
        terra = Group(eixo,earth)
        terra.rotate(OBLIQUIDADE*DEG,axis=UP)
        self.frame.reorient(0,80,0,earth.get_center(),CELESTIAL_SPHERE_RADIUS*3),
        self.wait(5)
        t= ValueTracker(0)
        earth.add_updater(lambda m:m.move_to([np.cos(t.get_value()/360*2*PI)* RENDER_SUN_EARTH_DISTANCE,np.sin(t.get_value()/360*2*PI)* RENDER_SUN_EARTH_DISTANCE,0]))
        earth.add_updater(lambda m,dt: m.rotate(dt*20,axis=(eixo.get_end()-eixo.get_start())))
        eixo.add_updater(lambda m:m.move_to(earth.get_center()))
        self.frame.add_updater(lambda m:m.reorient(0,80,0,earth.get_center(),CELESTIAL_SPHERE_RADIUS*3))
        self.play(t.animate.increment_value(360),run_time=20,rate_func=linear)
        earth.clear_updaters()
        eixo.clear_updaters()
        earth.add_updater(lambda m,dt: m.rotate(-dt*0.3,axis=OUT,about_point=earth.get_center()))
        eixo.add_updater(lambda m,dt: m.rotate(-dt*0.3,axis=OUT,about_point=earth.get_center()))
        earth.add_updater(lambda m,dt: m.rotate(dt*20,axis=(eixo.get_end()-eixo.get_start())))
        self.wait(10)
        

        

class Dedao(Scene):
    def construct(self):
        stars_data = extract_star_data()
        stars = Group()
        for x, y, z, size, color in stars_data:
            star = GlowDots(
                points=[2*np.array([x, y, z])],
                color=color,
                radius=2*size,
                opacity=0.75
            )
            stars.add(star)
        self.add(stars)

        terra = Earth(clouds=False)
        self.frame.reorient(0,80,0,terra.get_center(),CELESTIAL_SPHERE_RADIUS*1),
        self.add(terra)
        eixo = EixoPolar(90)
        self.play(ShowCreation(eixo))
        terra.add_updater(lambda m,dt: m.rotate(dt*0.2,axis=OUT))

class CircumpolarDedo(Scene):
    def construct(self):
        self.frame.reorient(0,90,0,ORIGIN,CELESTIAL_SPHERE_RADIUS*1.5)
        stars_data = extract_star_data(send_name=True)
        stars = Group()
        star_lookup = {}
        for star_dic in stars_data:
            x, y, z = star_dic["x"], star_dic["y"], star_dic["z"]
            size = star_dic["size"]
            color = star_dic["color"]
            hip = star_dic["hip"]
            distance = np.linalg.norm([x, y, z])
            
            star_lookup[hip] = star_dic
            star = GlowDots(
                points=[np.array([x, y, z])*60/distance],
                color=color,
                radius=size*200/distance,
                opacity=0.75
            )
            
            stars.add(star)
        self.add(stars)
        ASTERISM_LINES = [
            (11767, 77055),   
            (79822, 77055),    
            (72607, 77055),    
            (75097, 72607),   
            (75097, 79822),  
            (67301, 65378),  # Alkaid – Mizar
            (65378, 62956),  # Mizar – Alioth
            (62956, 59774),  # Alioth – Megrez
            (59774, 54061),  # Megrez – Dubhe
            (59774, 58001),  # Megrez – Phecda
            (58001, 53910),  # Phecda – Merak
            (53910, 54061),
            (60718, 61084),
            (62434, 59747),
            (107089,112405),
            (70638,107089),
            (70638,112405)
            
  


        ]

        # === 3. Criar as linhas dos asterismos ===
        lines = Group()
        for hip1, hip2 in ASTERISM_LINES:
            if hip1 not in star_lookup or hip2 not in star_lookup:
                continue

            s1 = star_lookup[hip1]
            s2 = star_lookup[hip2]

            p1 = np.array([s1["x"], s1["y"], s1["z"]]) * 60 / np.linalg.norm([s1["x"], s1["y"], s1["z"]])
            p2 = np.array([s2["x"], s2["y"], s2["z"]]) * 60 / np.linalg.norm([s2["x"], s2["y"], s2["z"]])

            line = Line3D(
                p1, p2,
                color=WHITE,
                opacity=0.7
            )
            lines.add(line)


        earth = Earth(clouds=False)
        self.add(earth)
        eixo = EixoPolar(90,comprimento=140,espessura=0.04*ELEMENTS_SCALE)

        earth.add_updater(lambda m,dt: m.rotate(0.1*dt,axis=OUT))
        self.wait(2)
        self.play(self.frame.animate.reorient(0,156,0,ORIGIN,CELESTIAL_SPHERE_RADIUS*1.5),run_time=2)
        self.play(ShowCreation(eixo))
        self.wait(2)
        self.play(self.frame.animate.add_ambient_rotation(0.1),run_time=2)
        self.frame.add_ambient_rotation(0.1)
        self.wait(2)
        earth.add_updater(lambda m:self.add(earth))
        self.play(earth.animate.set_opacity(0.4))
        self.play(self.frame.animate.reorient(phi_degrees=175,gamma_degrees=0,center=ORIGIN+[0,0,16],height=CELESTIAL_SPHERE_RADIUS*1.2),run_time=4)
        self.wait(6)
        self.play(self.frame.animate.reorient(phi_degrees=90,gamma_degrees=0,center=ORIGIN,height=CELESTIAL_SPHERE_RADIUS*8),run_time=4)
        self.wait(28)
        self.play(earth.animate.set_opacity(1))
 
class LuaEfeito(Scene):
    def construct(self):
        stars_data = extract_star_data()
        stars = Group()
        for x, y, z, size, color in stars_data:
            star = GlowDots(
                points=[np.array([x, y, z])],
                color=color,
                radius=2*size,
                opacity=0.75
            )
            stars.add(star)
        self.add(stars)
        
        self.frame.reorient(0,90,0,ORIGIN,CELESTIAL_SPHERE_RADIUS*1.8)
        earth = Earth(clouds=False)
        self.add(earth)
        moon = Moon()
        self.add(moon)
        moon.shift(RIGHT*CELESTIAL_SPHERE_RADIUS*1.5)
        earth.add_updater(lambda m,dt: m.rotate(0.7*dt,axis=OUT))
        moon.add_updater(lambda m,dt: m.rotate(0.06*dt,axis=OUT,about_point=ORIGIN))
        moon.add_updater(lambda m,dt: m.rotate(0.06*dt,axis=OUT))
        self.wait(10)
        self.play(self.frame.animate.reorient(0,0,0,ORIGIN,CELESTIAL_SPHERE_RADIUS*3.4),run_time=2)
        self.wait(6)
        self.play(self.frame.animate.add_ambient_rotation(0.7),run_time=4,rate_func=linear)
        self.frame.add_ambient_rotation(0.7)
        self.wait(15)

class TerraOrbitandoSol(Scene):
    def construct(self):
        # Fundo de estrelas
        stars_data = extract_star_data()
        stars = Group()
        for x, y, z, size, color in stars_data:
            star = GlowDots(
                points=[np.array([x, y, z])],
                color=color,
                radius=2 * size,
                opacity=0.75
            )
            stars.add(star)
        self.add(stars)

        # Reorientar a câmera inicial
        self.frame.reorient(0, 90, 0, ORIGIN, CELESTIAL_SPHERE_RADIUS * 1.8)

        # Criar o Sol
        sun = Sun()
        self.add(sun)

        # Criar a Terra
        earth = Earth(clouds=False).scale(0.6)
        self.add(earth)

        # Colocar a Terra na posição inicial da órbita
        earth.shift(RIGHT * CELESTIAL_SPHERE_RADIUS * 3)

        # Rotação própria da Terra
        earth.add_updater(lambda m, dt: m.rotate(0.7 * dt, axis=OUT))

        # Translação orbital da Terra ao redor do Sol
        earth.add_updater(lambda m, dt: m.rotate(0.07 * dt, axis=OUT, about_point=ORIGIN))

        # Esperar e mover a câmera suavemente
        self.wait(10)
        self.play(
            self.frame.animate.reorient(0, 0, 0, ORIGIN, CELESTIAL_SPHERE_RADIUS * 8),
            run_time=2
        )
        self.wait(10)
        self.play(
            self.frame.animate.reorient(0, 0, 0, [-6.25605774, 57.06007385 , 0.        ], CELESTIAL_SPHERE_RADIUS * 8),run_time=2
            )
        print(earth.get_center())
        self.frame.add_updater(lambda m: m.reorient(phi_degrees=0, gamma_degrees=0, center=earth.get_center(), height=CELESTIAL_SPHERE_RADIUS * 8))
        self.wait(10)
        self.frame.add_ambient_rotation(0.7)
        self.wait(15)



class PrecessaoFinal(Scene):
    def construct(self):
        stars_data = extract_star_data()
        stars = Group()
        for x, y, z, size, color in stars_data:
            star = GlowDots(
                points=[2*np.array([x, y, z])],
                color=color,
                radius=4*size,
                opacity=0.75
            )
            stars.add(star)

        
        self.add(stars)
        
        sun=Sun(radius=RENDER_SUN_RADIUS*5)
        mercury = Mercury().shift(RIGHT * RENDER_SUN_MERCURY_DISTANCE)
        venus = Venus().shift(RIGHT * RENDER_SUN_VENUS_DISTANCE)
        earth = Earth(clouds=False).shift(RIGHT * RENDER_SUN_EARTH_DISTANCE)
        mars = Mars().shift(RIGHT * RENDER_SUN_MARS_DISTANCE)
        jupiter = Jupiter().shift(RIGHT * RENDER_SUN_JUPITER_DISTANCE)
        saturn = Saturn().shift(RIGHT * RENDER_SUN_SATURN_DISTANCE)
        uranus = Uranus().shift(RIGHT * RENDER_SUN_URANUS_DISTANCE)
        neptune = Neptune().shift(RIGHT * RENDER_SUN_NEPTUNE_DISTANCE)

        # ========== OPCIONAL: ÓRBITAS ==========
        orbits = [
            Circle(radius=RENDER_SUN_MERCURY_DISTANCE, color=GREY_A, stroke_opacity=0.2),
            Circle(radius=RENDER_SUN_VENUS_DISTANCE, color=GREY_A, stroke_opacity=0.2),
            Circle(radius=RENDER_SUN_EARTH_DISTANCE, color=GREY_A, stroke_opacity=0.2),
            Circle(radius=RENDER_SUN_MARS_DISTANCE, color=GREY_A, stroke_opacity=0.2),
            Circle(radius=RENDER_SUN_JUPITER_DISTANCE, color=GREY_A, stroke_opacity=0.2),
            Circle(radius=RENDER_SUN_SATURN_DISTANCE, color=GREY_A, stroke_opacity=0.2),
            Circle(radius=RENDER_SUN_URANUS_DISTANCE, color=GREY_A, stroke_opacity=0.2),
            Circle(radius=RENDER_SUN_NEPTUNE_DISTANCE, color=GREY_A, stroke_opacity=0.2),
        ]
        sistema = Group()
        for orbit in orbits:
            sistema.add(orbit)
        # ========== ADICIONA OS CORPOS NA CENA ==========
        sistema.add(
            sun,
            mercury, venus, mars,
            jupiter, saturn, uranus, neptune
        )
        # sistema.shift(LEFT * RENDER_SUN_EARTH_DISTANCE)
        # earth.shift(LEFT * RENDER_SUN_EARTH_DISTANCE)
        self.add(sistema)
        self.add(earth)
        
        
        eixo=EixoPolar(90).move_to(earth.get_center())
        self.add(eixo)
        
        sun.move_to(ORIGIN)
        earth.add_updater(lambda m,dt:m.rotate(dt*0.1,axis=eixo.get_start()-eixo.get_end(),about_point=m.get_center()))
        terra = Group(eixo,earth)
        terra.rotate(OBLIQUIDADE*DEG,axis=UP)
        t= ValueTracker(-90)
        earth.move_to([np.cos(t.get_value()/360*2*PI)* RENDER_SUN_EARTH_DISTANCE,np.sin(t.get_value()/360*2*PI)* RENDER_SUN_EARTH_DISTANCE,0])
        eixo.add_updater(lambda m:m.move_to(earth.get_center()))
        self.frame.reorient(0,80,0,earth.get_center(),CELESTIAL_SPHERE_RADIUS*20),
        self.wait(5)
        
        
        earth.add_updater(lambda m:m.move_to([np.cos(t.get_value()/360*2*PI)* RENDER_SUN_EARTH_DISTANCE,np.sin(t.get_value()/360*2*PI)* RENDER_SUN_EARTH_DISTANCE,0]))
        earth.add_updater(lambda m,dt: m.rotate(dt*20,axis=(eixo.get_end()-eixo.get_start())))
        self.play(t.animate.increment_value(360),run_time=20,rate_func=linear)
        t2 = ValueTracker(1)
        self.frame.add_updater(lambda m:m.reorient(0,80,0,[np.cos(t.get_value()/360*2*PI)* RENDER_SUN_EARTH_DISTANCE,np.sin(t.get_value()/360*2*PI)* RENDER_SUN_EARTH_DISTANCE,0],CELESTIAL_SPHERE_RADIUS*20*t2.get_value()))
        self.play(t.animate.increment_value(180),t2.animate.increment_value(-0.9),rate_func=linear,run_time=20)
        earth.add_updater(lambda m,dt: m.rotate(-dt*0.3,axis=OUT,about_point=earth.get_center()))
        eixo.add_updater(lambda m,dt: m.rotate(-dt*0.3,axis=OUT,about_point=earth.get_center()))
        earth.add_updater(lambda m,dt: m.rotate(dt*20,axis=(eixo.get_end()-eixo.get_start())))
        self.play(t.animate.increment_value(100001),rate_func=linear,run_time=15)
        
        self.wait(10)
        
    


        
# Example usage
class VideoScene(Scene):
    def construct(self):

        video = VideoMobject("StarTrail", frame_rate=30, loop=True)
        video.set_height( FRAME_HEIGHT*1)

        self.add(video)
        self.wait(5)  # Let the video play for 10 seconds
        # celestia_sphere = EsferaCeleste()
        # celestia_sphere.add_updater(lambda m: self.add(celestia_sphere))
        # self.add(celestia_sphere)
        superficie = SuperficieObservador()
        Norte = Text("N",font_size=700).next_to(superficie.get_corner(UP),20* UP)
        Leste = Text("L",font_size=700).next_to(superficie.get_corner(RIGHT),20* RIGHT)
        Sul = Text("S",font_size=700).next_to(superficie.get_corner(DOWN),20* DOWN)
        Oeste = Text("O",font_size=700).next_to(superficie.get_corner(LEFT),20* LEFT)
        stars,_=Stars(sphere_mode=True,glow=False,sphere_radius=RENDER_CELESTIAL_SPHERE_RADIUS,size_factor=1)
        group = Group(superficie,Norte,Leste,Oeste,Sul)
        group.rotate(-90*DEG,axis=RIGHT).rotate(-180*DEG,axis=UP).shift(RENDER_CELESTIAL_SPHERE_RADIUS*OUT+0.1*RENDER_CELESTIAL_SPHERE_RADIUS*DOWN)
        stars.rotate((8.62)*DEG,axis=RIGHT).shift(superficie.get_center())
        self.add(stars)
        superficie.apply_depth_test()
        stars.apply_depth_test()
        stars.add_updater(lambda m,dt:m.rotate(dt*0.1,axis=[0,0.15,-1],about_point=RENDER_CELESTIAL_SPHERE_RADIUS*OUT+0.1*RENDER_CELESTIAL_SPHERE_RADIUS*DOWN))
        self.play(FadeIn(superficie),FadeIn(Norte),FadeIn(Leste),FadeIn(Sul),FadeIn(Oeste),FadeIn(stars))
        arrow = Line3D([0,0,0],20*np.array([0,0.15,-1])).shift(superficie.get_center())
        self.play(FadeIn(arrow))
        self.wait(5)
        self.play(self.frame.animate.reorient(0,-10,0,ORIGIN,RENDER_CELESTIAL_SPHERE_RADIUS*3))
        self.wait(8)

class OlhandoNorte(Scene):
    def construct(self):
        # Posiciona a câmera olhando para o Norte (eixo Y positivo)
        # com uma leve inclinação para uma boa visualização.
        self.frame.reorient(0, 85, 0, ORIGIN+[0,0,RENDER_CELESTIAL_SPHERE_RADIUS * 0.1], RENDER_CELESTIAL_SPHERE_RADIUS * 0.5)

        # Cria os elementos da cena usando suas classes
        esfera_celeste = EsferaCeleste()
    
        superficie = SuperficieObservador()
        stars,_=Stars(sphere_mode=True,glow=False,sphere_radius=RENDER_CELESTIAL_SPHERE_RADIUS,size_factor=1)
        stars.rotate((-90)*DEG,axis=RIGHT)
        Norte = Text("N",font_size=700).next_to(superficie.get_corner(UP),20* UP)
        Leste = Text("L",font_size=700).next_to(superficie.get_corner(RIGHT),20* RIGHT)
        Sul = Text("S",font_size=700).next_to(superficie.get_corner(DOWN),20* DOWN)
        Oeste = Text("O",font_size=700).next_to(superficie.get_corner(LEFT),20* LEFT)
        # Adiciona os elementos à cena
        self.play(
            ShowCreation(esfera_celeste),
            ShowCreation(superficie),
            Write(Norte),
            Write(Leste),
            Write(Sul),
            Write(Oeste)
        )
        esfera_celeste.add_updater(lambda m: self.add(esfera_celeste))
        
        self.wait(5)
        
        self.play(self.frame.animate.reorient(0, 70, 0, ORIGIN, RENDER_CELESTIAL_SPHERE_RADIUS * 2.5))

        self.play(FadeIn(stars))
        moon = Moon(radius=RENDER_CELESTIAL_SPHERE_RADIUS*0.05).shift(RIGHT*RENDER_CELESTIAL_SPHERE_RADIUS+OUT*RENDER_CELESTIAL_SPHERE_RADIUS*0.1)
        sun = Sun(radius=RENDER_CELESTIAL_SPHERE_RADIUS*0.05).move_to(IN*RENDER_CELESTIAL_SPHERE_RADIUS)
        self.play(ShowCreation(moon),ShowCreation(sun))
        self.wait(2)
        stars.add_updater(lambda m,dt:m.rotate(dt*3*-0.1,axis=[0,1,0],about_point=ORIGIN))
        moon.add_updater(lambda m,dt:m.rotate(dt*3*-0.06,axis=[0,1,0],about_point=ORIGIN))
        sun.add_updater(lambda m,dt:m.rotate(dt*3*-0.1,axis=[0,1,0],about_point=ORIGIN))
        self.wait(15)