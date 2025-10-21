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

class OLAA(Scene):
    def construct(self):
        self.frame.reorient(90,45,0,ORIGIN,CELESTIAL_SPHERE_RADIUS*2)
        light = self.camera.light_source
        esfera_celeste = EsferaCeleste()
        superficie = SuperficieObservador()
        Norte = Tex("N",font_size=7).next_to(superficie.get_corner(UP),0.2* UP)
        Leste = Tex("E",font_size=7).next_to(superficie.get_corner(RIGHT),0.2* RIGHT)
        Sul = Tex("S",font_size=7).next_to(superficie.get_corner(DOWN),0.2* DOWN)
        Oeste = Tex("O",font_size=7).next_to(superficie.get_corner(LEFT),0.2* LEFT)
        self.add(Norte,Leste,Sul,Oeste)
        self.add(esfera_celeste, superficie)
        esfera_celeste.add_updater(lambda m: self.add(esfera_celeste))
        # self.embed()
        
        eixo = EixoPolar(-22.9)
        self.play(ShowCreation(eixo))
        
        
        marcador = MarcadorAnguloExterno(PontoAstro(22.9,180),PontoAstro(0,180),cor=PURPLE,espessura=4)
        latitude = Tex(r'\varphi',font_size=3).rotate(90*DEGREES,axis=RIGHT).rotate(90*DEGREES,axis=OUT).move_to(marcador.get_center()*1.1)
        self.play(ShowCreation(marcador),Write(latitude))
        
        
        equador = Equador(-22.9)
        self.play(ShowCreation(equador))
        
        meridianolocal = MeridianoLocal()
        self.play(ShowCreation(meridianolocal))
        
        sol = Sun(radius=0.008,texture="data_image/sun.jpg").move_to(P(51,0).get_center())
        self.play(FadeIn(sol))
        polo = PontoAstroEquatorial(66.56,270,-22.9,43)
        tsl = ValueTracker(43)
        ecliptica = GrandeCirculo(polo.get_center())
        self.play(ShowCreation(ecliptica))
        
        ecliptica.add_updater(lambda m:m.become(GrandeCirculo(PontoAstroEquatorial(66.56,270,-22.9,tsl.get_value()).get_center())))
        self.play(tsl.animate.increment_value(94))
        self.play(tsl.animate.increment_value(-94))
        ecliptica.clear_updaters()
        self.embed()

class Gurjao(Scene):
    def construct(self):
        #USAR True em Vmobject
        self.frame.reorient(90,45,0,ORIGIN,CELESTIAL_SPHERE_RADIUS*2)
        light = self.camera.light_source
        esfera_celeste = EsferaCeleste()
        superficie = SuperficieObservador()
        Norte = Tex("N",font_size=7).next_to(superficie.get_corner(UP),0.2* UP)
        Leste = Tex("E",font_size=7).next_to(superficie.get_corner(RIGHT),0.2* RIGHT)
        Sul = Tex("S",font_size=7).next_to(superficie.get_corner(DOWN),0.2* DOWN)
        Oeste = Tex("O",font_size=7).next_to(superficie.get_corner(LEFT),0.2* LEFT)
        self.add(Norte,Leste,Sul,Oeste)
        self.add(esfera_celeste, superficie)
        esfera_celeste.add_updater(lambda m: self.add(esfera_celeste))
        
        ponto1 = P(0,90,cor="#FF00FF")
        lugar_geometrico_1 = GrandeCirculo(ponto1.get_center())
        seta = Seta(ponto1,0,cor="#FF00FF",comprimento=0.04,espessura=0.2)
        angle = ValueTracker(0)
        seta.add_updater(lambda m: m.become(Seta(ponto1,angle.get_value(),cor="#FF00FF",comprimento=0.04,espessura=0.2)))
        
        #Latituded
        self.save_state()

        eixo = EixoPolar(42)
        self.play(ShowCreation(eixo))
        
        #Ponto
        self.save_state()
    
        self.play(ShowCreation(ponto1))
        

        
        
        #Qualquer direação 
        self.add(seta)
        marcador =  MarcadorAngulo(ponto1,lugar_geometrico_1.point_from_proportion(angle.get_value()/360),raio_arco=0.05)
        marcador.add_updater(lambda m: m.become(MarcadorAngulo(ponto1,lugar_geometrico_1.point_from_proportion(angle.get_value()/360),raio_arco=0.05)))
        self.play(angle.animate.increment_value(360),run_time=5,rate_func=linear)
        
        #Marcador
        self.add(marcador)
        
        #equador
        self.play(ShowCreation(lugar_geometrico_1),angle.animate.increment_value(360),run_time=8,rate_func=linear)
        self.remove(marcador,seta)
        ponto2 = P(0,90,cor="#FF00FF")
        self.add(ponto2)
        
        #sun pra RA
        sun = Sun(radius=0.07*ELEMENTS_SCALE, big_glow_ratio=0.1*RENDER_SUN_RADIUS).move_to(PontoAstroEquatorial(0,111.3,42,-2.5*15+111.3).get_center())
        sun2 = Sun(radius=0.07*ELEMENTS_SCALE, big_glow_ratio=0.1*RENDER_SUN_RADIUS).move_to(PontoAstroEquatorial(0,0,42,0).get_center())
        self.play(ShowCreation(sun))
        ar = Seta(sun,240,cor="#2C00F3",comprimento=0.04,espessura=0.2)
        self.play(ShowCreation(ar))
        
        #ANgulo horario do sol 
        angulo_sol = GrandeArco(sun,sun2)
        self.play(ShowCreation(angulo_sol))
        
        #conta
        eq2 = Tex(r"\alpha_{sol} - H_{sol} = 7:25 - 1:28 = 5:57").shift(DOWN*3).deactivate_depth_test().scale(0.8)

        eq2.fix_in_frame()
        eq2.set_z_index(6)
        
        self.play(Write(eq2))
        #7:25 a ascenção reta do sol - H(d0 sol de 10:32 -> 1:28)  então da a ascenção reta de 5:57 
        
        #tirar
        self.remove(angulo_sol,sun,eq2,ar)        
        
        
        
        #Outra observação
        self.play(ponto2.animate.move_to(P(0,130).get_center()))
        lugar_geometrico_2 = GrandeCirculo(ponto2.get_center(),cor=GREEN_SCREEN)
        marcador_40 = MarcadorAngulo(ponto2,ponto1,raio_arco=0.05)
        self.play(ShowCreation(marcador_40))

        #Segundo grande arco
        self.play(ShowCreation(lugar_geometrico_2))
        
        #equador e constelaa
        estrelas = Group()
        
        sco_stars = [[-26.43194444, 247.3519583, 0.91], 
             [-37.10375, 263.4022083, 1.62], 
             [-42.99783333, 264.3297083, 1.86], 
             [-22.62161111, 240.083375, 2.29], 
             [-34.29261111, 252.5426667, 2.29],
             [-39.02991667, 265.622, 2.39], 
             [-19.80538889, 241.3592917, 2.62],
             [-37.29575, 262.691, 2.7],
             [-28.21597222, 248.9706667, 2.82], 
             [-26.11405556, 239.713, 2.89], 
             [-25.59275, 245.2971667, 2.9], 
             [-40.12697222, 266.8961667, 2.99],
             [-38.04733333, 252.9676667, 3.0], 
             [-37.04336111, 267.464375, 3.19], 
             [-43.2385, 258.03825, 3.32], 
             [-38.01747222, 253.0839583, 3.56],
             [-42.36075, 253.6462917, 3.62],
             [-29.214, 239.2212083, 3.87], 
             [-20.66913889, 241.7017917, 3.93], 
             [-19.46063889, 242.9989167, 4.0], 
             [-35.25536111, 249.0935833, 4.18], 
             [-34.70433333, 247.8455833, 4.24], 
             [-38.63486111, 264.136875, 4.26], 
             [-20.86866667, 241.85125, 4.31], 
             [-24.16927778, 245.1590833, 4.55], 
             [-27.92630556, 243.075875, 4.58], 
             [-25.32708333, 238.4030417, 4.59], 
             [-25.75122222, 237.7447917, 4.63], 
             [-42.36202778, 253.498875, 4.7]]

        cru_stars = [[-59.68, 191.9305, 1.25],

            [-63.080094444, 186.64975, 1.4],

            [-57.1, 187.791375, 1.6],

            [-63.830055556, 186.6520833, 2.09],

            [-58.73111111, 183.7865, 2.79],

            [-60.44863889, 185.340875, 3.59],

            ]


        for c in range(0,18):
            star = sco_stars[c]
            estrelas.add(PontoAstroEquatorial(star[0], star[1],latitude=0,TSL_graus=200, tamanho=DEFAULT_DOT_RADIUS/12*np.e**(-0.33*star[2])))
        for star in cru_stars:
            estrelas.add(PontoAstroEquatorial(star[0], star[1],latitude=0,TSL_graus=200, tamanho=DEFAULT_DOT_RADIUS/12*np.e**(-0.33*star[2])))

        equador = Equador(42,cor=BLUE)
        self.play(ShowCreation(equador), ShowCreation(estrelas))
        
        #passado o tempo
        self.remove(marcador_40)
        self.play(estrelas.animate.rotate(30*DEGREES,axis=eixo.get_end()-eixo.get_start(),about_point=ORIGIN),lugar_geometrico_1.animate.rotate(30*DEGREES,axis=eixo.get_end()-eixo.get_start()),ponto1.animate.rotate(30*DEGREES,axis=eixo.get_end()-eixo.get_start(),about_point=ORIGIN),run_time=4)
        
        #tirar equador
        self.play(FadeOut(equador))
        self.play(FadeOut(estrelas))
        
        #inter
        inter = PontoAstroEquatorial(-1.263656537,89.25,42,119.25)
        self.play(ShowCreation(inter))
        
        #p
        zenite = P(90,0)
        polo = P(42,0)
        self.play(FadeIn(zenite),FadeIn(polo))
        
        arco1 = GrandeArco(zenite,polo,espessura=LINE_SIZE*3,cor=BLUE_E,num_pontos=200).deactivate_depth_test()
        arco2 = GrandeArco(inter,polo,espessura=LINE_SIZE*3,cor=BLUE_E,num_pontos=200).deactivate_depth_test()
        arco3 = GrandeArco(inter,zenite,espessura=LINE_SIZE*3,cor=BLUE_E,num_pontos=200).deactivate_depth_test()
        self.play(ShowCreation(arco1),ShowCreation(arco2),ShowCreation(arco3))
        angulo_interno1 = AnguloEsferico(polo,zenite,inter,raio_circulo=0.03,reducing_factor=60)
        angulo_interno2 = AnguloEsferico(zenite,inter,polo,raio_circulo=0.015)
        self.play(ShowCreation(angulo_interno1),ShowCreation(angulo_interno2))
        
        
        
        #labels
        arco_1_label = Tex(r'90^{\circ}-\varphi',font_size=1)
        arco_1_label.move_to(arco1.get_center()+OUT*0.025+RIGHT*0.008)
        arco_1_label.set_z_index(4)
        arco_1_label.rotate(-80*DEGREES,axis=OUT)
        
        arco_2_label = Tex(r'90^{\circ}-\delta',font_size=1)
        arco_2_label.move_to(arco2.get_center_of_mass()+OUT*0.015+LEFT*0.03)
        arco_2_label.set_z_index(4)
        arco_2_label.rotate(-110*DEGREES,axis=OUT)
        
        angle_1_label=Tex(r'140^{\circ}',font_size=1)
        angle_1_label.move_to(angulo_interno2.get_center()+LEFT*0.015+UP*0.005)
        angle_1_label.set_z_index(4)
        angle_1_label.rotate(-100*DEGREES,axis=OUT)
        
        angle_2_label=Tex(r't',font_size=1)
        angle_2_label.move_to(angulo_interno1.get_center()+DOWN*0.015+UP*0.005+OUT*0.005+LEFT*0.005)
        angle_2_label.set_z_index(4)
        angle_2_label.rotate(-90*DEGREES,axis=OUT)
        
        #valores arcos
        self.play(Write(arco_1_label),Write(arco_2_label))
        #azimute
        self.play(Write(angle_1_label))
        #t
        self.play(Write(angle_2_label))
        
        eq1 = Tex(
            r"\cos ", r"(", r"90^{\circ}", r"-", r"\phi", r")", 
            r"\cos  ", r"( ", r"t", r") ", 
            r"=", 
            r"\sin ", r"(  ", r"90^{\circ}", r" - ", r"\phi", r")  ", 
            r"\cot ", r"(   ", r"90^{\circ}", r" -  ", r"\delta", r")   ", 
            r"-", 
            r"\sin  ", r"(    ", r"t", r")    ", 
            r"\cot  ", r"(     ", r"140^{\circ}", r")     "
        ).shift(DOWN*3).deactivate_depth_test().scale(0.8)

        eq1.fix_in_frame()
        eq1.set_z_index(6)
        
        self.play(Write(eq1))

class HeitorCarta(Scene):
    def construct(self):
        self.frame.reorient(90,45,0,ORIGIN,CELESTIAL_SPHERE_RADIUS*2)
        light = self.camera.light_source
        esfera_celeste = EsferaCeleste()
        superficie = SuperficieObservador()
        Norte = Tex("N",font_size=7).next_to(superficie.get_corner(UP),0.1* UP)
        Leste = Tex("E",font_size=7).next_to(superficie.get_corner(RIGHT),0.1* RIGHT)
        Sul = Tex("S",font_size=7).next_to(superficie.get_corner(DOWN),0.1* DOWN)
        Oeste = Tex("W",font_size=7).next_to(superficie.get_corner(LEFT),0.1* LEFT)
        self.add(Norte,Leste,Sul,Oeste)
        self.add(esfera_celeste, superficie)
        esfera_celeste.add_updater(lambda m: self.add(esfera_celeste))
    
        latitude = ValueTracker(-70)
        eixo = EixoPolar(latitude.get_value())
        self.embed()
        #explicação sem animação 
        
        #latitude genérica no sul
        self.play(ShowCreation(eixo))
        
        #Considerando o sol numa posição aleatória
        sol = Sun(radius=0.008).move_to(P(50,40).get_center())
        self.play(FadeIn(sol))
        
        #Zenite&azimute
        zenite =P(90,0,tamanho=0.5*ASTRO_SIZE)
        self.play(ShowCreation(zenite))
        arco_azimute =GrandeArco(zenite,sol)
        self.play(ShowCreation(arco_azimute))
        meridiano = MeridianoLocal()
        self.play(ShowCreation(meridiano))
        angulo_azimute = AnguloEsferico(zenite,P(0,0),sol,cor=BLUE,espessura=4)
        self.play(ShowCreation(angulo_azimute))
        azimute_label = Tex(r'A',font_size=1.5)
        azimute_label.move_to(angulo_azimute.get_center_of_mass()+OUT*0.003+UP*0.012+RIGHT*0.008)
        azimute_label.set_z_index(4)
        self.play(Write(azimute_label))
        
        #polo
        polo = P(-latitude.get_value(),180,cor=RED)
        self.play(ShowCreation(polo))
        declina =GrandeArco(polo,sol)
        self.play(ShowCreation(declina))
        angulo_horario = AnguloEsferico(polo,zenite,sol,cor=PURPLE,espessura=4)
        self.play(ShowCreation(angulo_horario))
        
        #Angulo horário
        angulo_horario.add_updater(lambda m: m.become(AnguloEsferico(polo,zenite,sol,cor=PURPLE,espessura=4)))
        declina.add_updater(lambda m: m.become(GrandeArco(polo,sol)))
        eixo.add_updater(lambda m : m.become(EixoPolar(latitude.get_value())))
        arco_azimute.add_updater(lambda m: m.become(GrandeArco(zenite,sol)))
        angulo_azimute.add_updater(lambda m: m.become(AnguloEsferico(zenite,P(0,0),sol,cor=BLUE,espessura=4)))
        polo.add_updater(lambda m: m.become(P(-latitude.get_value(),180,cor=RED)))
        
        self.play(latitude.animate.set_value(-90),run_time=1)
        self.play(latitude.animate.set_value(-20),run_time=1)
        self.play(latitude.animate.set_value(-40),run_time=1)
        
        angulo_horario.clear_updaters()
        declina.clear_updaters()
        eixo.clear_updaters()
        arco_azimute.clear_updaters()
        angulo_azimute.clear_updaters()
        polo.clear_updaters()
        
        #explicação
        #qualitativa do pq a unica coisa que importa é a diferença do angulo horario com o azimute e dai o resto entra nisso
        
        #b) 
        arco1 = GrandeArco(polo,sol,espessura=LINE_SIZE*3,cor=BLUE_E,num_pontos=200).deactivate_depth_test()
        arco2 = GrandeArco(zenite,sol,espessura=LINE_SIZE*3,cor=BLUE_E,num_pontos=200).deactivate_depth_test()
        arco3= GrandeArco(zenite,polo,espessura=LINE_SIZE*3,cor=BLUE_E,num_pontos=200).deactivate_depth_test()
        self.play(ShowCreation(arco1),ShowCreation(arco2),ShowCreation(arco3))
        
        angulo_horario_label = Tex(r'H',font_size=1.5)
        angulo_horario_label.move_to(angulo_horario.get_center_of_mass()+OUT*0.005+UP*0.016+RIGHT*0.004)
        angulo_horario_label.set_z_index(4)
        self.play(Write(angulo_horario_label))
        
        angulo_sup = AnguloEsferico(zenite,sol,polo,espessura=4,cor=YELLOW)
        self.play(ShowCreation(angulo_sup))
        
        sup_label = Tex(r'180^{\circ}-A',font_size=1.2)
        sup_label.move_to(angulo_sup.get_center_of_mass()+OUT*0.003+DOWN*0.012+RIGHT*0.01)
        sup_label.rotate(70*DEGREES,axis=OUT)
        sup_label.set_z_index(4)
        self.play(Write(sup_label))
        
        dec_label = Tex(r'90^{\circ}+\delta',font_size=1.3)
        dec_label.move_to(arco1.get_center_of_mass()+OUT*0.006+DOWN*0.012+RIGHT*0.008)
        dec_label.rotate(70*DEGREES,axis=OUT)
        dec_label.set_z_index(4)
        self.play(Write(dec_label))
        
        phi_label = Tex(r'90^{\circ}-|\varphi|',font_size=1.3)
        phi_label.move_to(arco3.get_center_of_mass()+OUT*0.006+DOWN*0.012+LEFT*0.008)
        phi_label.rotate(100*DEGREES,axis=OUT)
        phi_label.set_z_index(4)
        self.play(Write(phi_label))
        
        eq1= Tex(r"\cot(90^\circ + \delta)\cdot\sin(90^\circ - |\varphi|) = \cot(180^\circ - A)\cdot\sin(-H) + \cos(90^\circ - |\varphi|)\cdot\cos(-H)")
        eq1.fix_in_frame()
        eq1.set_z_index(6)
        eq1.scale(0.8).shift(DOWN*3.5)
        self.play(Write(eq1))
        
        eq2 = Tex(r"-",r"H",r"= ",r"A")
        eq3 = Tex(r"\varepsilon ",r"= ",r"H",r"+",r"A")
        eq4 = Tex(r"\varepsilon ",r"= ",r"H",r"+",r"\tan^{-1}\left( \frac{-\sin(H)}{\tan(\delta)\cdot\cos(|\varphi|) + \sin(|\varphi|)\cdot\cos(H)} \right)")
        
        eq2.fix_in_frame()
        eq2.set_z_index(6)
        eq2.scale(0.8).shift(UP*3.5)
        self.play(Write(eq2))
        
        eq3.fix_in_frame()
        eq3.set_z_index(6)
        eq3.scale(0.8).shift(UP*3.5)
        
        eq4.fix_in_frame()
        eq4.set_z_index(6)
        eq4.scale(0.8).shift(UP*3.5)
        
        self.wait(2)
        self.play(TransformMatchingStrings(eq2,eq3))
        self.wait(3)
        self.play(TransformMatchingStrings(eq3,eq4))
        self.wait(2)
        
class TimeEquation(Scene):
    def construct(self):
        self.frame.reorient(90,45,0,ORIGIN,CELESTIAL_SPHERE_RADIUS*2)
        terra = Earth(radius=2*RENDER_EARTH_RADIUS,day_texture="data_image/earth_day_high.jpg",night_texture="data_image/earth_night.jpg",clouds=False)
        esfera_celeste = EsferaCeleste()
        light = self.camera.light_source
        
        
        esfera_celeste.add_updater(lambda m: self.add(esfera_celeste))
        self.add(terra,esfera_celeste)
        hora = ValueTracker(11)
        hora_verdadeira = ValueTracker(13)
        self.frame.reorient(0, 90, 0, ORIGIN, CELESTIAL_SPHERE_RADIUS*2)
        
        observador = PontoAstroEquatorial(0,-72.89,90,raio=2*RENDER_EARTH_RADIUS,cor="#250B48")
        observador1 = PontoAstroEquatorial(-89,-72.89,90,cor="#250B48")
        observador2 = PontoAstroEquatorial(89,-72.89,90,cor="#250B48")
        observador3 = Sun(radius=0.07*ELEMENTS_SCALE, big_glow_ratio=0.1*RENDER_SUN_RADIUS).move_to(observador.get_center()/RENDER_EARTH_RADIUS/2*CELESTIAL_SPHERE_RADIUS)
        observador4 = Sun(radius=0.07*ELEMENTS_SCALE, big_glow_ratio=0.1*RENDER_SUN_RADIUS).move_to(-observador.get_center()/RENDER_EARTH_RADIUS/2*CELESTIAL_SPHERE_RADIUS)
        meridiano = GrandeArco(observador1,observador2,cor=YELLOW)
        sun = Sun(radius=0.07*ELEMENTS_SCALE, big_glow_ratio=0.1*RENDER_SUN_RADIUS).move_to(np.array([np.cos(hora.get_value()*15*DEGREES+PI/2-72.89*DEGREES),np.sin(hora.get_value()*15*DEGREES+PI/2-72.89*DEGREES),0])*CELESTIAL_SPHERE_RADIUS)
        sun_verdadeira = Sun(radius=0.07*ELEMENTS_SCALE, big_glow_ratio=0.1*RENDER_SUN_RADIUS).move_to(np.array([np.cos(hora_verdadeira.get_value()*15*DEGREES+PI/2-72.89*DEGREES),np.sin(hora_verdadeira.get_value()*15*DEGREES+PI/2-72.89*DEGREES),0])*CELESTIAL_SPHERE_RADIUS)
        sun_verdadeira.set_opacity(0.1)
        light.move_to(np.array([np.cos(hora.get_value()*15*DEGREES+PI/2-72.89*DEGREES),np.sin(hora.get_value()*15*DEGREES+PI/2-72.89*DEGREES),0])*CELESTIAL_SPHERE_RADIUS*10)
        self.add(sun,sun_verdadeira,meridiano,observador)
        meridiano_rodado = MeridianoLocal().rotate(+PI/2,axis=UP,about_point=ORIGIN).rotate(-PI/2+(90-72.89)*DEGREES,axis=OUT,about_point=ORIGIN)
        group = Group(observador,terra,meridiano,observador3,meridiano_rodado)
        
        angulo_example = GrandeArco(observador3,sun_verdadeira)
        self.add(angulo_example,meridiano_rodado)
        self.embed()
        self.remove(angulo_example)
        sun.move_to(observador3)
        self.play(group.animate.rotate(30*DEGREES,axis=OUT,about_point=ORIGIN),run_time=5)
        self.play(self.frame.animate.reorient(270,-72.89,0,1.2*observador.get_center(),CELESTIAL_SPHERE_RADIUS*0.1))
        self.frame.reorient(0, 90, 0, ORIGIN, CELESTIAL_SPHERE_RADIUS*2)

class Lua(Scene):
    def construct(self):
        self.frame.reorient(90,45,0,ORIGIN,0.1)

        self.frame.add_updater(lambda m,dt:m.increment_theta(dt*10*DEGREES))
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
        moon = Moon(radius=RENDER_MOON_RADIUS*4)
        self.add(moon)
        moon.move_to(P(0,0,raio=RENDER_CELESTIAL_SPHERE_RADIUS*9).get_center())
        moon.add_updater(lambda m,dt: m.rotate(dt*10*DEGREES,axis=OUT,about_point=ORIGIN))
        self.embed()
        moon.clear_updaters()
        moon.add_updater(lambda m,dt: m.rotate(dt*-10*DEGREES,axis=OUT,about_point=ORIGIN))
        self.wait(20)
        